pkg load control
close all
clear all
s = tf('s');
%2.2
%Visto o octave não simplificar expressões iremos usar a expressão simplificada a que chegamos
%na alínea 2.1
%Para K_theta pertencente a valores reais positivos
K_theta = 1;
OLF = K_theta/(s^2)
CLF = K_theta/(s^2 + K_theta)
%Temos o Root Locus
figure(1);
rlocus(OLF);
%Step response da closed loop function
figure(2);
step(CLF,100);
%Os polos encontram-se no eixo dos imaginários e à medida que aumentamos o valor de K_theta os polos afastam-se de zero


%2.4
%Visto o octave não simplificar expressões iremos usar a expressão simplificada a que chegamos
%na alínea 2.3
%Para K_w pertencente a valores reais positivos
K_w = 1;
OLF = K_w/(s)
CLF = K_w/(s + K_w)
%Temos o Root Locus
figure(3);
rlocus(OLF);
%Step response da closed loop function
figure(4);
step(CLF,100);
%O polo encontra-se no eixo real negativo e á medida que aumentamos o valor de K_theta o polo afasta-se de zero

%Visto o octave não simplificar expressões iremos usar a expressão simplificada a que chegamos
%na alínea 2.3
%Para K_theta pertencente a valores reais positivos
K_w = 1;
Ls_23 = K_w/(s);

%Vamos adicionar um polo arbitrário ao controlador tal que s = z_1 seja menor que zero
z_1 = 11;
Ks = (K_w*s)/(s + z_1); %novo controlador
%Pelos cálculos já realizados ficamos com o seguinte Open loop
OLF = K_w/(s*(s+z_1))
%Root locus da Open Loop
figure (3);
rlocus(OLF);
%E ficamos com a seguinte Closed Loop
CLF = K_w/(s*(s+z_1)+K_w)
%Step response da closed loop function
figure(4);
step(CLF,100);


