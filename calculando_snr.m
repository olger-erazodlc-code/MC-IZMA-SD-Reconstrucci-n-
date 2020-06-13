clear, close all
L=2^15;
t=1:L;

a=sin(t)+.5*cos(2*t)+2;

clear, close all
amplitud=2;
T=linspace(0,1,24)';
seno1=amplitud*sin(2*pi*9*T);
potencia_senal=var(seno1)+mean(seno1)^2

sigma_v=4;
ruido=sigma_v.*randn(size(seno1));
potencia_ruido=var(ruido)+mean(ruido)^2
if (potencia_ruido<potencia_senal)
db=potencia_senal/potencia_ruido
else
db=-(potencia_ruido/potencia_senal)
end