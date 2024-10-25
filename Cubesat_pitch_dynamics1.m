pkg load control
close all
clear all
b_z = 22960.2*10^-9; %Campo Magnético que passa pelo cubo
J_y =1.667*10^-3; %Momento de inércia segundo o eixo dos yy;
a = b_z/J_y; %Constante
s = tf('s');
G = a/(s^2); %Função transferência a analisar

%1.3
figure (1);
bode(G); %Diagrama de Bode da função tranferência
grid on;

%1.5
figure (2);
step(G,100); %Step response da função transferência
