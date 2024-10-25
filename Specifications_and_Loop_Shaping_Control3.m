pkg load control
close all
clear all
s = tf('s')
%3.3 e 3.4
%Primeira parte
Ti=10.5;
Tw=2.5;

Gs = 1/s^2
Ks = (1 + Tw*s + 1/(Ti*s));
OLP = Gs*Ks
figure(1);
rlocusx(OLP);

K_theta = 5;
CLP = minreal((K_theta*(OLP))/(K_theta*(OLP) + 1))
figure(2);
step(CLP);

figure(3);
margin(CLP);
