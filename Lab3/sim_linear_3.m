close all;
clear;
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

f1 = [0;sin(deg2rad(30));-cos(deg2rad(30))];
f2 = [0;sin(deg2rad(30));-cos(deg2rad(30))];
f3 = [0;sin(deg2rad(-30));-cos(deg2rad(-30))];
f4 = [0;sin(deg2rad(-30));-cos(deg2rad(-30))];

n1 = cross(p1,f1);
n2 = cross(p2,f2);
n3 = cross(p3,f3);
n4 = cross(p4,f4);

%Parametros de simulação
Dt = 0.01;
t = 0:Dt:60;
Nsim = length(t);
beta = -0.5*Cd*ro_mars*l_L*l_L;

%Linearização-equilibrio
px = 0;
py = 0;
pz = -21;
p = [px;py;pz];
wx = 0;  
wy = 0;
wz = 0;
omega = [wx;wy;wz]; 
vx = 0;
vy = 0;
vz = 0;
v = [vx;vy;vz];
phi = deg2rad(0);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];
T = (0.25*mt*g_mars)/cos(deg2rad(30));
u = ones(Nsim,1)*[T,T,T,T];

a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), [zeros(2,3);0 0 2*beta*vz/mt], skew(g_mars*z_I), skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/mt*f1, 1/mt*f2, 1/mt*f3, 1/mt*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);

figure(1);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'b');
legend('p_x','p_y','p_z','v_z')
grid on;

[V,DL,W] = eig(A); % W'*A = D*W'

%Teste de estabilidade, controlabilidade e observabilidade
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('O Sistema não é controlável'); end
if rank(obsv(A,C)) < size(A,1), disp('O Sistema não é observável'); end
ctrb_modes = W'*B,
obsv_modes = C*V,

%SVD
G = tf(sys);
%plot dos valores singulares
figure(2);
sigmaplot(G);
grid on;
title('Singular Values');

%LQR
Q = blkdiag(1,1,10,1,1,100,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
X = A-B*Klqr;
CLlqr = ss(X,B,C,D);
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

%Simulação do LQR
Dt = 0.01;
t = 0:Dt:100;
NSim = length(t);
r = [0;0;-21;0;0;0;0;0;0;0;0;0]*(t>=0);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-21;0;0;0;0;0;0;0;0;0];
%C = eye(12);
for k = 1:NSim

    y(:,k) = C*x(:,k);

    u(:,k) = -Klqr*[y(:,k) - r(:,k)];

    x_dot = A*x(:,k) + B*u(:,k); % system derivatives
    xp = x(:,k) + Dt*x_dot; % integrate system state
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(3);
subplot(3,1,1);
plot(t,r(2,:),t,r(3,:),t,y(1,:),t,y(2,:),t,y(3,:));
ylabel('Posição [m]');
xlabel('Tempo [s]');
legend('rp_y','rp_z','p_x','p_y','p_z');
ylim([-25 10]);
grid on;

subplot(3,1,2);
plot(t,r(6,:),t,y(4,:),t,y(5,:),t,y(6,:));
ylabel('Velocidade [m/s]');
xlabel('Tempo [s]');
legend('r_v_z','v_x','v_y','v_z');
grid on;

subplot(3,1,3);
plot(t,u(1,:),t,u(2,:),t,u(3,:),t,u(4,:));
ylabel('Thrust [N]');
xlabel('Tempo [s]');
legend('T_1','T_2','T_3','T_4');
grid on;

%H_inf
B1 = zeros(12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12);
D12 = [zeros(8,4);W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

A0 = A;
B0 = [B1,B2];
C0 = [C1;C2];
D0 = [D11 D12; D21 D22];

P = ss(A0,B0,C0,D0);

nmeas = 12;
ncont = 4;

% tests on P:
nx = size(A0,1);
if (nx - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (nx - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
if (size(D12,2) - rank(D12)) > 0, disp('A2.1 on P: D12 is column rank deficient'); else disp('A2.1 on P: OK'); end
if (size(D21,1) - rank(D21)) > 0, disp('A2.2 on P: D21 is row rank deficient'); else disp('A2.1 on P: OK'); end
syms w real; 
Aux1 = [A0 - j*w*eye(size(A0)) , B2 ; C1 , D12];
if (size(Aux1,2) - rank(Aux1)) > 0,  disp('A3 on P: matrix is column rank deficient'); else disp('A3 on P: OK'); end
Aux2 = [A0 - j*w*eye(size(A0)) , B1 ; C2 , D21];
if (size(Aux2,1) - rank(Aux2)) > 0,  disp('A4 on P: matrix is column rank deficient'); else disp('A4 on P: OK'); end

[Kinf,CLinf,gammainf,info_inf] = hinfsyn(P,nmeas,ncont);
poles_CLinf = pole(CLinf);
if any(real(poles_CLinf) >= 0), disp('CL system with Hinf controller not stable'); else, disp('CL system with Hinf controller is stable'); end

for k = 1:NSim

    y(:,k) = C*x(:,k);

    v = -[(r(:,k)-y(:,k))];
    u(:,k) = Kinf.C*v;

    x_dot = A*x(:,k) + B*u(:,k); 
    xp = x(:,k) + Dt*x_dot; 
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(4);
subplot(3,1,1);
plot(t,r(2,:),t,r(3,:),t,y(1,:),t,y(2,:),t,y(3,:));
ylabel('Posição [m]');
xlabel('Tempo [s]');
legend('rp_y','rp_z','p_x','p_y','p_z');
ylim([-25 10]);
grid on;

subplot(3,1,2);
plot(t,r(6,:),t,y(4,:),t,y(5,:),t,y(6,:));
ylabel('Velocidade [m/s]');
xlabel('Tempo [s]');
legend('r_v_z','v_x','v_y','v_z');
grid on;

subplot(3,1,3);
plot(t,u(1,:),t,u(2,:),t,u(3,:),t,u(4,:));
ylabel('Thrust [N]');
xlabel('Tempo [s]');
legend('T_1','T_2','T_3','T_4');
grid on;

figure(5);
sigmaplot(CLlqr,'k',Kinf);
grid on;
title('Sigma Values - CL (no weights)');
legend('LQR','H\infty');
