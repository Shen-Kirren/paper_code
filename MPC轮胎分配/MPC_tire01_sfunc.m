function [sys,x0,str,ts]=MPC_tire01_sfunc(t,x,u,flag)
% u�����������������MPCģ����x(k)����������[gama,beta,MFx,MFy,efl,efr,erl,err]��Ŀǰ��û���MFx/MFy��ô����
switch flag
    case 0
        % ��ʼ��
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 2
        % ״̬����
        sys=mdlUpdates(t,x,u);
    case 3
        % �������
        sys=mdlOutputs(t,x,u);
    case {1,4,9}
        % ��Ч����
        sys=[];
    otherwise
        % �������룬��flag��ӡ����
        error(['unhandled flag=',num2str(flag)]);
end
end

function [sys,x0,str,ts] = mdlInitializeSizes
% mdlInitializeSizes - ��ʼ������ϵͳ
%
% Syntax: [sys,x0,str,ts] = mdlInitializeSizes
%
% ��Ҫ���ǽ�x0�����ֵ���涨S-function����������ĸ��������͡��Լ�����/��ɢϵͳ����������
    sizes=simsizes;
    sizes.NumContStates=0;  % ����״̬���ĸ���
    sizes.NumDiscStates=8;  % ��ɢ״̬���ĸ���
    % sizes.NumContStates+sizes.NumDiscStates������������x��ά��
    sizes.NumOutputs=4;     % S-funtion������ĸ���
    sizes.NumInputs=14;      % S-funtion�������ĸ�����u��ά��,�����˺�ڽ��ٶ�gama,���Ĳ�ƫ��beta,�ĸ����ӵ�ת�٣����ٶ�ax,ay,�ĸ����ӵ�������
    sizes.DirFeedthrough=1; % ֱ�Ӵ��ݾ��󣬼�y=Cx+Du��D����
    sizes.NumSampleTimes=1;% ����ʱ��ĸ�������һ��������ʱ��
    sys=simsizes(sizes);   % �������ֵ��sys���
    x0=[0;0;0;0;0;0;0;0];             % ״̬��x�ĳ�ʼֵx0�����ó�ʼ��������ͺ�ڽǶ���0
    global U;               % ȫ�ֱ���U�����U����S-function��u����ϵͳ״̬���̵�U��������ϵͳX_dot=AX+BU��U
    U=[0.5;0.5;0.5;0.5];                % U�а������ǳ��ٺ�ת�ǣ�ά��Ϊ2
    str = [];               % Set str to an empty matrix.strΪ����������mathworks��˾��û�����ô������һ���ڳ�ʼ���н����Ϳռ���
    ts  = [0.1 0];          % tsΪ1*2ά�ȵ���������������+ƫ����sample time: [period, offset]�����濪ʼ0s��ÿ0.1������һ��
end

function sys = mdlUpdates(t,x,u)
%mdlUpdates - ״̬����
%
% Syntax: sys = mdlUpdates(t,x,u)
%
% ״̬����ʱ����ϵͳ��Ϊ��x_dot=Ax+Bu�е�x��Ϊ�������
    sys=x;
end

function sys = mdlOutputs(t,x,u)
%mdlOutputs - ���������������MPC����Ҫ���֣�������״̬���̵Ĺ��������ݺ����Ĺ������Լ�QP��������
%
% Syntax: sys = mdlOutputs(t,x,u)
%
% �����t,x,u��ǰ���Сдt,x,u����һ�£���S-function�Ĳ���
global u_piao r A_c B_c C_c E_c A_d B_d C_d E_d C
global kesi U
tic
Nx=8;                   %״̬����������ʵ��X����X~��ά�ȴ�С
Nu=4;                   %��������������ʵ��U����U~��ά�ȴ�С
Nw=6;
Ny=10;                  %�۲���������������gama`,beta`,ei`,Ti`
Np=6;                   %Ԥ��ʱ��6���������ڣ������в�������ΪTs=0.1s������predict horizonΪ0.6s
Nc=3;                   %����ʱ��3���������ڣ������в�������ΪTs=0.1s������control horizonΪ0.3s
Row=10;                 %�ɳ�����Ȩ��

% ���������Լ������������
T=0.1;
L=2.6;
lf=1.04;
lr=1.56;
bf=1.481;
br=1.486;
m=1111;
% Row=10;                 %�ɳ�����Ȩ��
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
re=0.304;               %��Ч�����뾶304mm
Ps1_factor=Ps1/max(U(1),0.5);
Ps2_factor=Ps2/max(U(2),0.5);
Ps3_factor=Ps3/max(U(3),0.5);
Ps4_factor=Ps4/max(U(4),0.5);

% �����µ�״̬��kesi
kesi=zeros(Nx+Nu,1);
Ks=m/2/(lf+lr)^2*(lr/Kf-lf/Kr); %Kf/Kr��ǰ���ֵĲ�ƫ�Ǹն�
gama_ref=min(vx/((lf+lr)+(1+Ks*vx^2)));
beta_ref=0;             %�ο������Ĳ�ƫ��Ϊ0
MFx_ref=0;
MFy_ref=0;
efl_ref=0; efr_ref=0; erl_ref=0; err_ref=0; %���ֵ�ת�ٲ�ֵ�Ĳο�Ϊ0����ϣ�����ֲ�Ҫ������
kesi(1)=u(1)-gama_ref;
kesi(2)=u(2)-beta_ref;
kesi(3)=0;
kesi(4)=0;
kesi(5)=u(5); kesi(6)=u(6); kesi(7)=u(7); kesi(8)=u(8);
kesi(9)=U(1); kesi(10)=U(2); kesi(11)=U(3); kesi(12)=U(4); % �൱�ڸ��ֵ�ת�زο�����0

u_piao=zeros(Nu,1);

% ����kesi(k+1)=A`*kesi(k)+B`*u_piao���µ�״̬�ռ���ʽ
% �ȹ���X(k+1)=A_d*X(k)+B_d*U(k)+E_d*W(k);
A_c=zeros(Nx);
% B_c=zeros(Nx,Nu);
A_c(1,3)=1/Iz;
A_c(1,4)=1/Iz;
A_c(2,1)=-1;
tao=0.1;               % tao:delay factor����ϵͳ�Ĺ���ʱ�ӣ�����ȡ0.1����Ϊ0.1����һ��
A_c(3,3)=1/tao;
lf_factor=1+((vy+lf*gama)/vx)^2;
lr_factor=1+((vy-lf*gama)/vx)^2;
factor=lf/vx;
alpha_dot=[-factor/lf_factor, -factor/lf_factor, factor/lr_factor, factor/lr_factor]'.*gama_dot;
ky=diag(Cyfl, Cyfr, Cyrl, Cyrr); % Cy����̥��ƫ�Ǹն�,Fy=Cy*alpha
Ay=[bf/2*sin(delta_f)+lf*cos(delta_f), -bf/2*sin(delta_f)+lf*cos(delta_f), -lr, -lr];
km=Ay*ky*alpha_dot;
A_c(4,3)=km/Iz;
A_c(4,4)=km/Iz;
re_factor=2*Iz*re;      % re����Ч�����뾶
A_c(5,3)=-bf/re_factor;
A_c(5,4)=-bf/re_factor;
A_c(6,3)=bf/re_factor;
A_c(6,4)=bf/re_factor;
A_c(7,3)=-br/re_factor;
A_c(7,4)=-br/re_factor;
A_c(8,3)=br/re_factor;
A_c(8,4)=br/re_factor;

re_factor=2*Iz*re;      % re����Ч�����뾶
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


% ���潫���������״̬��״̬�ռ���ʽ��ɢ������
% ��õ���״̬���Ʊ��ʽ
% X(k+1)=A_d(k)*X(k)+B_d(k)*U(k)+E_d(k)*W(k);
% Y(k+1)=C_d(k+1)*X(k+1);
% ��Ҫ˵����A/B/C/D_d���󲢷�ʱ����ģ��������Ų�ͬʱ�̲�����仯�ģ���Ϊ�������һЩ�������������lf,lr���־�bf,br֮���
% ������һЩ�仯���������糵��vx,vy,���ٶ�ax,ay,��ڽ��ٶ�gama,���Ĳ�ƫ��beta
dt=0.05;                % ��������/�������0.05s��20Hz
A_d=eye(Nx)+A_c.*dt;
B_d=B_c.*dt;
C_d=C_c;
E_d=E_c.*dt;

% ����A`,B`
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
% %ȥ����ei�Ĺ۲⣬��Ϊei�ǵ�i�����Ӵ�ʱ������һʱ��ת�ٵĲ�ֵ��û����
C(3,9)=1*Ps1_factor; C(4,10)=1*Ps2_factor; C(5,11)=1*Ps3_factor; C(6,12)=1*Ps4_factor;

% �µ�ϵͳԤ��������̣�Y=PHI*kesi+THETA*U��p86ʽ4.12
% ����PHI��THETA����Ҳ��ʹ��cell2mat�ķ���
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

Q=100*eye(Nx*Np,Nx*Np); %״̬����Ȩ�ؾ���
R=50*eye(Nu*Nc);        %��������Ȩ�ؾ���

% Ŀ�꺯��J=y'*Py*y+u_piao'*Pu*u_piao
% �����u_piao���Ǧ�u`��Ҳ����MPC����Ż�������X
% ��Ϊy=PHI*kesi+THETA*u_piao
% ��Jչ����������E=PHI*kesi,�ã�
% J=E'*Py*E+2*E'*Py*THETA*u_piao+u_piao'*(THETA'*Py*THETA+Pu)*u_piao
% ��u_piao��ΪMPC�е�X���������ɳ����ӣ�X��ά��Ϊu_piao+1
% ��ΪMPC���QP����ĸ�ʽ��minJ=1/2*X'*H*X+f'*X
H_cell{1,1}=THETA'*Q*THETA+R;%H�����(1,1)Ϊʲô��������ӣ���Ҫȥ��p85ʽ4.7�������p86ʽ4.10��4.12
H_cell{1,2}=zeros(Nu*Nc,1);
H_cell{2,1}=zeros(1,Nu*Nc);
H_cell{2,2}=Row;        %�ɳ����ӵ�Ȩ��ϵ��
H=cell2mat(H_cell);
H=(H+H')/2;

E=PHI*kesi;
f_cell=cell(1,2); %�˴�����ĸ�ʽ�Ѿ���f'����ʽ�ˣ����Բ���Ҫת��
f_cell{1,1}=2*E'*Q*THETA;
f_cell{1,2}=0;
f=cell2mat(f_cell);

% ����Լ������
% Լ��������������ת�ص����½磬ÿ�����������½磬�ɳ����ӵ����½�
du_min=-10; du_max=10; % ����ת�ص���������
u_min=0; u_max=150;    % ����ת�ص�����
% U_last=U;              % ά��Ϊ4x1����������һʱ�̵Ŀ������

A_t=tril(ones(Nc));
A_I=kron(A_t,eye(Nu));
Ut=kron(ones(Nc,1),U); %ÿ��ʱ�̵�U��Ҫ�ܵ����ƣ���������kron
Umin=kron(ones(Nc,1),ones(Nu,1).*u_min);
Umax=kron(ones(Nc,1),ones(Nu,1).*u_max);

% A_cons�а����������֣��ֱ���A_I*X<=Umax-Ut
% �Լ�-A_I*X<=-Umin+Ut,��A_I*X>=Umax+Ut
% �������cell�����ľ���ά��������(Nc*Nu)�ķ����ұ�һ��������������ƴ�����ģ�������Ϊ�˶�Ӧ�Ͼ���ά�ȣ���ΪX��ά�����ɳ����ӣ���Nc*Nu+1ά��������
A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1)};
b_cons_cell={Umax-Ut;-Umin+Ut};%�˴�Ut����δ���£�����Ut����һʱ�̵Ŀ�������������ĺ����Ǹ�ʱ�̵Ŀ������������ֵ��ΧӦ������һʱ�̿����������һ�����²���������
A_cons=cell2mat(A_cons_cell);
b_cons=cell2mat(b_cons_cell);
% ��Ϊ����MPC����������X������ʵ�ʿ���ϵͳ�п�������������������ʵ�ʿ���ϵͳ�еĿ�����������һʱ�̿���������һʱ�̿�����������ģ�Լ������������b_cons�Ĺ�����
M=10;                   %�ɳ���������
delta_Umin=kron(ones(Nc,1),ones(Nu,1).*du_min);
delta_Umax=kron(ones(Nc,1),ones(Nu,1).*du_max);
lb=[delta_Umin;0];
ub=[delta_Umax;M];

% ���QP���⣬ʹ�õ����ڵ㷨
options = optimset('Algorithm','interior-point-convex');
% �����ٴ�˵��һ�£�X���2x1�������������ʵ������������������ģ���Ƶ������е�u_piao��Ҳ����U-Uref
% ��Ϊ���ǲ�û������U�г��ٺ�ת�ǵķ��������������ʼ��U����0
[X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
% ��U���£�U+=u_piao,����Ҫ���ѵ�һ�����kesi(9:12,1)������Ԫ�ؾ���U
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