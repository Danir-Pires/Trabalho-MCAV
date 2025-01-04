%% Controlability and Observability
%% pkg load control

m = 1539;
zI=[0;0;1]

A = [
  zeros(3), eye(3), [0, -89, 0; 89, 0, 0; 0, 0, 0], zeros(3);
  zeros(3), -1.3268e-4 * ones(3), [0, -3.71, 0; 3.71, 0, 0; 0, 0, 0], [0, -89, 0; 89, 0, 0; 0, 0, 0];
  zeros(3), zeros(3), zeros(3), eye(3);
  zeros(3), zeros(3), zeros(3), zeros(3)
];

B = [
    zeros(3, 4);
    (1/m) * [
        [0; 0.342; -0.940], [0; 0.342; -0.940], [0; -0.342; -0.940], [0; -0.342; -0.940]
    ];
    zeros(3, 4);
    [
        [6.63e-4; -4.97e-4; 2.21e-4], [6.63e-4; 4.97e-4; 2.21e-4], [-6.63e-4; 4.97e-4; -2.21e-4], [-6.63e-4; -4.97e-4; 2.21e-4]
    ]
];
C = [   eye(3)    , zeros(3)   , zeros(3) , zeros(3)
        zeros(1,3), zeros(1,3) , zI'      , zeros(1,3)  ];
D = zeros(4);

sys = ss(A,B,C,D);

[Vj,Jor] = schur(A),
% [Vj,Jor] = jordan(sym(A)),
[V,DL,W] = eig(A),
% mode_obs = C*V
% mode_ctrl = W'*B

control_resultMatrix = zeros(1, 12);
obs_resultMatrix = zeros(1, 12);

disp('Controlability - mode analysis')
%checking if system is controllable
for i = 1:12
    aux=W(:,i)'*B; %grau de contrabilidade dos modos (só é preciso que uma entrada seja controlavel) - interpretar quais as entradas que são
    disp(i)
    disp(aux)
end

disp('Observability - mode analysis')
%checking if system is observable
for i = 1:12
    observable=true;
    obs_resultMatrix = C*V(:,i);
    disp(i)
    disp(obs_resultMatrix)
end

if any(real(diag(DL)) >=0 ), disp('Linearized system 1 is not stable.'); end
n1_unstable_modes = rank(ctrb(A,B))-12;
if n1_unstable_modes > 0, disp('Linearized system 1 is not controlable.'); end
n1_unobservable_modes = rank(obsv(A,C))-12;
if n1_unobservable_modes > 0, disp('Linearized system 1 is not observable.'); end

% compute transfer functionq
G = tf(ss(A,B,C,0));

% get LQR controller
Q = blkdiag(9*eye(3),10*eye(3),0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

%---%
A0 = A;
B1 = zeros(12,12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12,12);
D12 = [zeros(8,4); W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

B0 = [B1, B2];
C0 = [C1; C2];
D0 = [D11, D12; D21, D22];

P = ss(A0, B0, C0, D0);

nmeas = 12; % Number of measured outputs
ncont = 4; % Number of control inputs

% tests on P:
if (12 - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (12 - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
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

% Hinf controller
Dt = 0.01;
t = 0:Dt:300;
r = [0;0;-500;0;0;45;0;0;0;0;0;0]*(t>=0); %onde queremos que ele vá (referencia) OP1b
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-2100;0;0;89;0;0;0;0;0;0]; %onde ele está condicoes iniciais OP1a
C=eye(12)
for k = 1:NSim

    % get measurements:
    y1(:,k) = C*x(:,k);

    % get control action:
    
    v = -[(r(:,k)-y1(:,k))];
    u(:,k) = Kinf.C*v; % approximation for LQR-like performance of Hinf
    
    
    % simulate system:
    x_dot = A*x(:,k) + B*u(:,k); % system derivatives
    xp = x(:,k) + Dt*x_dot; % integrate system state
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(9041002);
plot(t,u);
grid on;
xlabel('$$t [s]$$');
ylabel('$$u(t)$$');
legend('$$u_1$$','$$u_2$$');

figure(9041003);
plot(t,y1(3,:),t,r(3,:)); %referencia em z - 3
grid on;
xlabel('$$t [s]$$');
legend('$$y_1$$','$$r_1$$');

