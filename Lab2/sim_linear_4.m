close all;

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
T = ((mt*g_mars*0.25)/cos(deg2rad(30)));
u = ones(Nsim,1)*[T,T,T,T];
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
vy = 25;
vz = 0;
v = [vx;vy;vz];
phi = deg2rad(10);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];

a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

c = [0, -g_mars, 0
     g_mars*cos(phi), 0, 0
     -g_mars*sin(phi), 0, 0];


A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), zeros(3), c, skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/mt*f1, 1/mt*f2, 1/mt*f3, 1/mt*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

C = [eye(3), zeros(3), zeros(3), zeros(3)
     zeros(3), eye(3), zeros(3), zeros(3)];

D = zeros(6,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);

figure(1);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'y',t,y_L(:,5),'b',t,y_L(:,6),'c');
legend('px','py','pz','vx','vy','vz')
grid on;

[V,D,W] = eig(A); % W'*A = D*W'

% Teste de estabilidade, controlabilidade e observabilidade
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('O Sistema não é controlável'); end
if rank(obsv(A,C)) < size(A,1), disp('O Sistema não é observável'); end
ctrb_modes = W'*B,
obsv_modes = C*V,
