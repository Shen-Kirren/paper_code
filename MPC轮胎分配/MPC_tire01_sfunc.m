function [sys,x0,str,ts]=MPC_tire01_sfunc(t,x,u,flag)
% u中输入的向量，就是MPC模型中x(k)的量，包括[gama,beta,MFx,MFy,efl,efr,erl,err]。目前还没想好MFx/MFy怎么计算
switch flag
    case 0
        % 初始化
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 2
        % 状态更新
        sys=mdlUpdates(t,x,u);
    case 3
        % 计算输出
        sys=mdlOutputs(t,x,u);
    case {1,4,9}
        % 无效操作
        sys=[];
    otherwise
        % 错误输入，把flag打印出来
        error(['unhandled flag=',num2str(flag)]);
end
end

function [sys,x0,str,ts] = mdlInitializeSizes
% mdlInitializeSizes - 初始化整个系统
%
% Syntax: [sys,x0,str,ts] = mdlInitializeSizes
%
% 主要就是将x0赋予初值，规定S-function的输入输出的个数、类型、以及连续/离散系统的其他设置
    sizes=simsizes;
    sizes.NumContStates=0;  % 连续状态量的个数
    sizes.NumDiscStates=8;  % 离散状态量的个数
    % sizes.NumContStates+sizes.NumDiscStates就是主函数中x的维度
    sizes.NumOutputs=4;     % S-funtion输出量的个数
    sizes.NumInputs=14;      % S-funtion输入量的个数，u的维度,包括了横摆角速度gama,质心侧偏角beta,四个轮子的转速，加速度ax,ay,四个轮子的纵向力
    sizes.DirFeedthrough=1; % 直接传递矩阵，即y=Cx+Du的D矩阵
    sizes.NumSampleTimes=1;% 采样时间的个数，就一个主采样时间
    sys=simsizes(sizes);   % 设置完后赋值给sys输出
    x0=[0;0;0;0;0;0;0;0];             % 状态量x的初始值x0，设置初始横纵坐标和横摆角都是0
    global U;               % 全局变量U，这个U不是S-function的u，是系统状态方程的U，即整个系统X_dot=AX+BU的U
    U=[0.5;0.5;0.5;0.5];                % U中包含的是车速和转角，维度为2
    str = [];               % Set str to an empty matrix.str为保留参数，mathworks公司还没想好怎么用它，一般在初始化中将其滞空即可
    ts  = [0.1 0];          % ts为1*2维度的向量，采样周期+偏移量sample time: [period, offset]，仿真开始0s后每0.1秒运行一次
end

function sys = mdlUpdates(t,x,u)
%mdlUpdates - 状态更新
%
% Syntax: sys = mdlUpdates(t,x,u)
%
% 状态更新时，将系统认为的x_dot=Ax+Bu中的x作为输出给出
    sys=x;
end

function sys = mdlOutputs(t,x,u)
%mdlOutputs - 计算输出，这里是MPC的主要部分，包括了状态方程的构建，传递函数的构建，以及QP问题的求解
%
% Syntax: sys = mdlOutputs(t,x,u)
%
% 这里的t,x,u和前面的小写t,x,u保持一致，是S-function的参量
global u_piao r A_c B_c C_c E_c A_d B_d C_d E_d C
global kesi U
tic
Nx=8;                   %状态量个数，其实是X或者X~的维度大小
Nu=4;                   %控制量个数，其实是U或者U~的维度大小
Nw=6;
Ny=10;                  %观测量个数，包括，gama`,beta`,ei`,Ti`
Np=6;                   %预测时域，6个采样周期，本例中采样周期为Ts=0.1s，所以predict horizon为0.6s
Nc=3;                   %控制时域，3个采样周期，本例中采样周期为Ts=0.1s，所以control horizon为0.3s
Row=10;                 %松弛因子权重

% 车辆参数以及仿真参数设置
T=0.1;
L=2.6;
lf=1.04;
lr=1.56;
bf=1.481;
br=1.486;
m=1111;
% Row=10;                 %松弛因子权重
Iz=2031.4;
Kf=370.667*2;
Kr=370.667*2;
Cyfl=370.667; Cyfr=370.667; Cyrl=370.667; Cyrr=370.667;
gama=u(1);
beta=u(2);
vx=u(3);
vy=u(4);
delta_f=u(9);
Ps1=u(10);
Ps2=u(11);
Ps3=u(12);
Ps4=u(13);
gama_dot=u(14);
re=0.304;               %有效滚动半径304mm
Ps1_factor=Ps1/max(U(1),0.5);
Ps2_factor=Ps2/max(U(2),0.5);
Ps3_factor=Ps3/max(U(3),0.5);
Ps4_factor=Ps4/max(U(4),0.5);

% 构造新的状态量kesi
kesi=zeros(Nx+Nu,1);
Ks=m/2/(lf+lr)^2*(lr/Kf-lf/Kr); %Kf/Kr是前后轮的侧偏角刚度
gama_ref=min(vx/((lf+lr)+(1+Ks*vx^2)));
beta_ref=0;             %参考的质心侧偏角为0
MFx_ref=0;
MFy_ref=0;
efl_ref=0; efr_ref=0; erl_ref=0; err_ref=0; %各轮的转速差值的参考为0，即希望各轮不要发生打滑
kesi(1)=u(1)-gama_ref;
kesi(2)=u(2)-beta_ref;
kesi(3)=0;
kesi(4)=0;
kesi(5)=u(5); kesi(6)=u(6); kesi(7)=u(7); kesi(8)=u(8);
kesi(9)=U(1); kesi(10)=U(2); kesi(11)=U(3); kesi(12)=U(4); % 相当于各轮的转矩参考量是0

u_piao=zeros(Nu,1);

% 构造kesi(k+1)=A`*kesi(k)+B`*u_piao的新的状态空间表达式
% 先构造X(k+1)=A_d*X(k)+B_d*U(k)+E_d*W(k);
A_c=zeros(Nx);
% B_c=zeros(Nx,Nu);
A_c(1,3)=1/Iz;
A_c(1,4)=1/Iz;
A_c(2,1)=-1;
tao=0.1;               % tao:delay factor控制系统的惯性时延，这里取0.1是因为0.1采样一次
A_c(3,3)=1/tao;
lf_factor=1+((vy+lf*gama)/vx)^2;
lr_factor=1+((vy-lf*gama)/vx)^2;
factor=lf/vx;
alpha_dot=[-factor/lf_factor, -factor/lf_factor, factor/lr_factor, factor/lr_factor]'.*gama_dot;
ky=diag(Cyfl, Cyfr, Cyrl, Cyrr); % Cy是轮胎侧偏角刚度,Fy=Cy*alpha
Ay=[bf/2*sin(delta_f)+lf*cos(delta_f), -bf/2*sin(delta_f)+lf*cos(delta_f), -lr, -lr];
km=Ay*ky*alpha_dot;
A_c(4,3)=km/Iz;
A_c(4,4)=km/Iz;
re_factor=2*Iz*re;      % re是有效滚动半径
A_c(5,3)=-bf/re_factor;
A_c(5,4)=-bf/re_factor;
A_c(6,3)=bf/re_factor;
A_c(6,4)=bf/re_factor;
A_c(7,3)=-br/re_factor;
A_c(7,4)=-br/re_factor;
A_c(8,3)=br/re_factor;
A_c(8,4)=br/re_factor;

re_factor=2*Iz*re;      % re是有效滚动半径
B_c=[-bf/re_factor, bf/re_factor, -br/re_factor, br/re_factor;
    0, 0, 0, 0;
    0, 0, 0, 0;
    -km*bf/re_factor, km*bf/re_factor, -km*br/re_factor, km*br/re_factor,;
    -1/Iw+bf*bf/(2*re*re_factor), -bf*bf/(2*re*re_factor), bf*br/(2*re*re_factor), -bf*br/(2*re*re_factor);
    -bf*bf/(2*re*re_factor), -1/Iw+bf*bf/(2*re*re_factor), -bf*br/(2*re*re_factor), bf*br/(2*re*re_factor);
    bf*br/(2*re*re_factor), -bf*br/(2*re*re_factor), -1/Iw+br*br/(2*re*re_factor), -br*br/(2*re*re_factor);
    -bf*br/(2*re*re_factor), bf*br/(2*re*re_factor), -br*br/(2*re*re_factor), -1/Iw+br*br/(2*re*re_factor);];
% C_c_cell=cell(2,2);
% C_c_cell{1,1}=[1, 0, 0, 0;
%                0, 1, 0, 0;];
% C_c_cell{1,2}=zeros(2,4);
% C_c_cell{2,1}=zeros(4);
% C_c_cell{2,2}=eye(4);
% C_c=cell2mat(C_c_cell);
E_c=[0, 0, 0, 0, 0, 0;
     0, 1/vx, 0, 0, 0, 0;
     0, 0, Ay(1)/tao, Ay(2)/tao, Ay(3)/tao, Ay(4)/tao;
     0, 0, 0, 0, 0, 0;
     1/re, 0, re/Iw, 0, 0, 0;
     1/re, 0, 0, re/Iw, 0, 0;
     1/re, 0, 0, 0, re/Iw, 0;
     1/re, 0, 0, 0, 0, re/Iw;];


% 下面将上面的连续状态的状态空间表达式离散化处理
% 获得的是状态递推表达式
% X(k+1)=A_d(k)*X(k)+B_d(k)*U(k)+E_d(k)*W(k);
% Y(k+1)=C_d(k+1)*X(k+1);
% 需要说明，A/B/C/D_d矩阵并非时不变的，而是随着不同时刻采样会变化的，因为里面除了一些不会变的例如轴距lf,lr，轮距bf,br之类的
% 还会有一些变化的量，比如车速vx,vy,加速度ax,ay,横摆角速度gama,质心侧偏角beta
dt=0.05;                % 采样周期/采样间隔0.05s，20Hz
A_d=eye(Nx)+A_c.*dt;
B_d=B_c.*dt;
C_d=C_c;
E_d=E_c.*dt;

% 构造A`,B`
A_cell=cell(2,2);
A_cell{1,1}=A_d;
A_cell{1,2}=B_d;
A_cell{2,1}=zeros(Nu,Nx);
A_cell{2,2}=eyes(Nu);
A=cell2mat(A_cell);

B_cell=cell(2,1);
B_cell{1,1}=B_d;
B_cell{2,1}=eye(Nu);
B=cell2mat(B_cell);

% E_cell=cell(2,1);
% E_cell{1,1}=E_d;
% E_cell{2,1}=zeros(Nx+Nu-Nw,Nw);
% E=cell2mat(E_cell);

C=zeros(6,12);
C(1,1)=1; C(2,2)=1; 
% C(3,5)=1; C(4,6)=1; C(5,7)=1; C(6,8)=1;
% %去掉对ei的观测，因为ei是第i个轮子此时刻与上一时刻转速的差值，没有用
C(3,9)=1*Ps1_factor; C(4,10)=1*Ps2_factor; C(5,11)=1*Ps3_factor; C(6,12)=1*Ps4_factor;

% 新的系统预测输出方程，Y=PHI*kesi+THETA*U，p86式4.12
% 构造PHI和THETA矩阵，也是使用cell2mat的方法
PHI_cell=cell(Np,1);
THETA_cell=cell(Np,Nc);
for j=1:1:Np
    PHI_cell{j,1}=C*A^j;
    for k=1:1:Nc
        if k<=j
            THETA_cell{j,k}=C*A^(j-k)*B;
        else
            THETA_cell{j,k}=zeros(Nx,Nu);
        end
    end
end
PHI=cell2mat(PHI_cell);
THETA=cell2mat(THETA_cell);

Q=100*eye(Nx*Np,Nx*Np); %状态量的权重矩阵
R=50*eye(Nu*Nc);        %控制量的权重矩阵

% 目标函数J=y'*Py*y+u_piao'*Pu*u_piao
% 这里的u_piao就是Δu`，也就是MPC最后优化出来的X
% 因为y=PHI*kesi+THETA*u_piao
% 将J展开，并且令E=PHI*kesi,得：
% J=E'*Py*E+2*E'*Py*THETA*u_piao+u_piao'*(THETA'*Py*THETA+Pu)*u_piao
% 将u_piao换为MPC中的X，考虑上松弛因子，X的维度为u_piao+1
% 改为MPC求解QP问题的格式，minJ=1/2*X'*H*X+f'*X
H_cell{1,1}=THETA'*Q*THETA+R;%H矩阵的(1,1)为什么是这个样子，需要去看p85式4.7，并结合p86式4.10与4.12
H_cell{1,2}=zeros(Nu*Nc,1);
H_cell{2,1}=zeros(1,Nu*Nc);
H_cell{2,2}=Row;        %松弛因子的权重系数
H=cell2mat(H_cell);
H=(H+H')/2;

E=PHI*kesi;
f_cell=cell(1,2); %此处构造的格式已经是f'的形式了，所以不需要转置
f_cell{1,1}=2*E'*Q*THETA;
f_cell{1,2}=0;
f=cell2mat(f_cell);

% 构造约束条件
% 约束条件包括各轮转矩的上下界，每次增量的上下界，松弛因子的上下界
du_min=-10; du_max=10; % 控制转矩的增量限制
u_min=0; u_max=150;    % 控制转矩的限制
% U_last=U;              % 维度为4x1，包括了上一时刻的控制输出

A_t=tril(ones(Nc));
A_I=kron(A_t,eye(Nu));
Ut=kron(ones(Nc,1),U); %每个时刻的U都要受到限制，所以用了kron
Umin=kron(ones(Nc,1),ones(Nu,1).*u_min);
Umax=kron(ones(Nc,1),ones(Nu,1).*u_max);

% A_cons中包含了两部分，分别是A_I*X<=Umax-Ut
% 以及-A_I*X<=-Umin+Ut,即A_I*X>=Umax+Ut
% 所以这个cell包含的矩阵维度是两个(Nc*Nu)的方阵，右边一列零向量，上下拼起来的，零向量为了对应上矩阵维度，因为X的维度有松弛因子，是Nc*Nu+1维的列向量
A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1)};
b_cons_cell={Umax-Ut;-Umin+Ut};%此处Ut矩阵还未更新，所以Ut是上一时刻的控制量矩阵，这里的含义是该时刻的控制量矩阵的数值范围应该在上一时刻控制量矩阵的一个上下波动区间内
A_cons=cell2mat(A_cons_cell);
b_cons=cell2mat(b_cons_cell);
% 因为我们MPC中求解出来的X是我们实际控制系统中控制量的增量，而我们实际控制系统中的控制量是用上一时刻控制量加这一时刻控制增量计算的，约束里体现在了b_cons的构造中
M=10;                   %松弛因子上限
delta_Umin=kron(ones(Nc,1),ones(Nu,1).*du_min);
delta_Umax=kron(ones(Nc,1),ones(Nu,1).*du_max);
lb=[delta_Umin;0];
ub=[delta_Umax;M];

% 求解QP问题，使用的是内点法
options = optimset('Algorithm','interior-point-convex');
% 这里再次说明一下，X这个2x1的列向量解出来实际上是我们整个控制模型推导过程中的u_piao，也就是U-Uref
% 因为我们并没有引入U中车速和转角的反馈，所以这个初始的U就是0
[X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
% 将U更新，U+=u_piao,还需要提醒的一点就是kesi(9:12,1)这两个元素就是U
u_piao(1)=X(1);
u_piao(2)=X(2);
u_piao(3)=X(3);
u_piao(4)=X(4);
U(1)=kesi(9)+u_piao(1);
U(2)=kesi(10)+u_piao(2);
U(3)=kesi(11)+u_piao(3);
U(4)=kesi(12)+u_piao(4);
u_real=U;
sys=u_real;
toc
end