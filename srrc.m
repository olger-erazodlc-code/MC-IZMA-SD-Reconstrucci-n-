function [response]=srrc(factor_sobrem,roll_off)
%factor_sobrem=factor sobremuestreo
%factor_sobrem=factor de sobremuestreo en este caso el numero de muestras
% es =10
% roll_off=0,5
a=roll_off; % en otros documentos B=roll_off=0,5
t=-4:1/factor_sobrem:4; %Limitando respuesta de -4T a 4T
% Este se incrementa o decrementa de acuerdo a los requerimientos
% t es el periodo, cuando t es grande existe la posibilidad de una
% menor interferencia intersombolica
%t=-4:1/factor_sobrem:4;  t=-4:0.1:4;  Tiene 80 elementos
%_____________________ cantidad de elementos en el periodo________________
 disp('Conjunto de elementos de t.'); 
 yy=length(t);
 disp(yy);
%________________ fin mostrarcantidad de elementos en el periodo___________
p=zeros(1,length(t));
    for i=1:1:length(t)
        if t(i)==0
            p(i)= (1-a)+4*a/pi; % falta multiplicar por 1/(sqrt(t))
        else if t(i)==1/(4*a) || t(i)==-1/(4*a)
               p(i)=a/sqrt(2)*((1+2/pi)*sin(pi/(4*a))+(1-2/pi)*cos(pi/(4*a)));
               % falto al segmento  a/sqrt(2) hacerlo a/sqrt(2*Ts) Ts=t
              else
                p(i) = (sin(pi*t(i)*(1-a))+4*a*t(i).*cos(pi*t(i)*(1+a)))./(pi*t(i).*(1-(4*a*t(i)).^2));
                % falta multiplicar por 1/sqrt(Ts)
             end
        end
    end
    % tal parece que Ts se asume como 1
    % ya que Ts es el recíproco de la velocidad de símbolos
    
    response=p./sqrt(sum(p.^2)); %Normalization to unit energy
    %response=matriz p  / raiz cuadrada de (sumatoria de (p al cuadrado) )
    % como algunos valores  son negativos, al elevarlos al cuadrado se
    % hacen positivos, despues se suman todos y se vuelve a sacar la raiz
    % cuadradd para normalizar p
    
end
