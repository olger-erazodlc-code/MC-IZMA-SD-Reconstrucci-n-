
function decibel = snr(x,y)
%snr - signal to noise ratio
%decibel = snr(x,y);
%decibel = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )

%   x es la señal original, referencia. Incluye ruido
%   y es la señal sin ruido
%  decibel = 20*log10(norm(x(:))/norm(x(:)-y(:)));
% 102 es señal con ruido/102 señal con ruido -100 sin ruido=102/(102-100)=100/2=50
% decibel =20*log10(50) =20*1.69=33

% ________________ LO SIGUIENTE SE IMPLEMENTA______________________
% Y si solo se compara la señal sin ruido, con  respecto al ruido?
decibel = 10*log10(norm(x(:))/norm(y(:)));
% x es la señal sin ruido, y es el ruido
% 10*log (100/2)=10*log10(50)=16,98
% 10*log (100/50)=10*log10(2)=3,01 db en este caso indica que la señal 
% esta al doble del ruido.
% 10*log (50/100)=10*log10(0.5)=-3,01 db en este caso indica que la señal 
% esta a la mitad de la señal del ruido.
end