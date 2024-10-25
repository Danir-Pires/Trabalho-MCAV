try %octave only
  pkg load control
end
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
rlocusx(OLF);
%Step response da closed loop function
figure(2);
step(CLF,100);
%Os polos encontram-se no eixo dos imaginários e à medida que aumentamos o valor de K_theta os polos afastam-se de zero


%2.4
%Visto o octave não simplificar expressões iremos usar a expressão simplificada a que chegamos
%na alínea 2.3
%Para K_theta pertencente a valores reais positivos
K_w = 1;
OLF = K_w/(s)
CLF = K_w/(s + K_w)
%Temos o Root Locus
figure(3);
rlocusx(OLF);
%Step response da closed loop function
figure(4);
step(CLF,100);
%O polo encontra-se no eixo real negativo e á medida que aumentamos o valor de K_theta o polo afasta-se de zero


%2.5
%Valores atrbuidos a K's
K_theta = 1; %Ganho que será atribuido a partir do Root Locus
     %Ganhos que irão variar o Root Locus e a resposta do Closed Loop
K_w = 7;
K_i = 0.013;

%Open Loop Function
OLF = (K_w*(s^2) + K_theta*s + K_i)/(s^3)
%CLosed Loop Function
CLF = (K_w*s^2 + K_theta*s + K_i)/(s^3 + K_w*s^2 + K_theta*s + K_i);

%Root Locus
figure(1);
rlocusx(OLF);

%Step Response
K_theta = 3; %Damos valores a K_theta e verificamos o step response
CLF = (K_w*s^2 + K_theta*s + K_i)/(s^3 + K_w*s^2 + K_theta*s + K_i)
figure(2);
step(CLF);


