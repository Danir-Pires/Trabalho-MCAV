clear 
close all
clc

% Parâmetros
g_mars = 3.721;
Cd = 1.05;
ro_mars = 0.02;
z_I = [0;0;1];

%Para o Lander
l_L = 3.0;
w_L = 2.7;
h_L = 2.2;
m = 1025; %Massa do lander é igual à massa do Rover
mt = 2*m;
J_L = [(1/12)*m*(l_L^2+h_L^2) 0 0  %Momento de Inércia do Lander
      0 (1/12)*m*(w_L^2+h_L^2) 0   
      0 0 (1/12)*m*(l_L^2+w_L^2)];

%Para o Rover
l_R = 2.7;
w_R = 3.0;
h_R = 2.2;

J_R = [(1/12)*m*(l_R^2+h_R^2) 0 0  %Momento de Inércia do Rover
      0 (1/12)*m*(w_R^2+h_R^2) 0
      0 0 (1/12)*m*(l_R^2+w_R^2)];

J = [J_L(1,1)+m*(1.1)^2+J_R(1,1)+m*(1.1)^2 0 0  %Momento de Inércia do Sistema
     0 J_L(2,2)+m*(1.1)^2+J_R(2,2)+m*(1.1)^2 0
     0 0 J_L(3,3)+J_R(3,3)];

%Posições dos retrorockets
p1 = [-1.35;-1.5;0];
p2 = [1.35;-1.5;0];
p3 = [1.35;1.5;0];
p4 = [-1.35;1.5;0];

%Thrust de cada retrorocket
%Como T1=T2=T3=T4 podemos admitir um T genérico tal que Tt=4*T
T = m*g_mars*0.25/cos(deg2rad(30));

%Vetor da Força Propulsiva
f1 = [0;sin(deg2rad(30));-cos(deg2rad(30))];
f2 = [0;sin(deg2rad(30));-cos(deg2rad(30))];
f3 = [0;sin(deg2rad(-30));-cos(deg2rad(-30))];
f4 = [0;sin(deg2rad(-30));-cos(deg2rad(-30))];
fp = T*f1+T*f2+T*f3+T*f4;


%Vetor de Momentos 
n1 = cross(p1,f1);
n2 = cross(p2,f2);
n3 = cross(p3,f3);
n4 = cross(p4,f4);
np = T*n1+T*n2+T*n3+T*n4;

% Simulação do Sistema não Linear
%Parâmetros
Dt = 0.01;
t = 0:Dt:60; %visto a fase de powered descend demorar aproximadamente 1 minuto
u = [fp;np]*ones(size(t)) - 0.5*[fp;np]*(t>=35);
nx = 12;
ny = 4;
x0 = [0;0;-2100;0;0;89;0;0;0;0;0;0];
Nsim = length(t);
%x = zeros(nx,Nsim);
%y = zeros(ny,Nsim);
x(:,1) = x0;

%Simulação
for k = 1:Nsim

    p = x(1:3,k);
    v = x(4:6,k);
    lbd = x(7:9,k);
    omega = x(10:12,k);
    R = Euler2R(lbd);
    Q = Euler2Q(lbd);
    fg = g_mars*R'*z_I;
    fa = -0.5*ro_mars*l_L*l_L*z_I;
    fp = u (1:3,k);
    np = u (4:6,k);

    p_dot = R*v;
    v_dot = -skew(omega)*v+1/mt*(fg+fa*v(3)^2+fp);
    lbd_dot = Q*omega;
    omega_dot = (-J^-1*skew(omega)*J*omega)+(J^-1*np);
    x_dot = [p_dot;v_dot;lbd_dot;omega_dot];

    x(:,k+1) = x(:,k) + Dt*x_dot;
  
end
x(:,k+1) = [];

folder = '..\figures\';
example_name = 'ex4-7_rocket3d_sim';
printfigs = 0;
vehicle3d_show_data(t,x,u,0,[folder example_name],90270);