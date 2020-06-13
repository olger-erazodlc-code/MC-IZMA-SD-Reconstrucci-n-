
function decibel = snr(x,y)
%snr - signal to noise ratio
%decibel = snr(x,y);
%decibel = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )

%   x es la se�al original, referencia. Incluye ruido
%   y es la se�al sin ruido
%  decibel = 20*log10(norm(x(:))/norm(x(:)-y(:)));
% 102 es se�al con ruido/102 se�al con ruido -100 sin ruido=102/(102-100)=100/2=50
% decibel =20*log10(50) =20*1.69=33

% ________________ LO SIGUIENTE SE IMPLEMENTA______________________
% Y si solo se compara la se�al sin ruido, con  respecto al ruido?
decibel = 10*log10(norm(x(:))/norm(y(:)));
% x es la se�al sin ruido, y es el ruido
% 10*log (100/2)=10*log10(50)=16,98
% 10*log (100/50)=10*log10(2)=3,01 db en este caso indica que la se�al 
% esta al doble del ruido.
% 10*log (50/100)=10*log10(0.5)=-3,01 db en este caso indica que la se�al 
% esta a la mitad de la se�al del ruido.
end