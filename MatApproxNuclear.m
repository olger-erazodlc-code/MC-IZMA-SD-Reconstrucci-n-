function [NewMAT, Error,cal_mat2] = MatApproxNuclear(A, Init, Zone, lambda, N, tol,con_desv, display)
% This function performs matrix approximation when including only some
% entries. In other words, it solves the problem NewMAT = argmin
% ||(NewMAT-A).*Zone||_F, s.t. ||NewMAT||_*<=lambda. Because of convexity
% the solution is global/optimal.
%
% Coded by Gil Shabat, 2012
% A=matriz a completar
% Init=matriz que se completa en cada llamado a función
% Zone=matriz de 0´s y 1´s, 0 son datos perdidos
% lambda 
% N es el número de iteraciones 
% tolerancia norma frobenius tol= 1e-7
% display=1, ya que si la norma es nuclear o la norma es espectral
% 1 muestra en pantalla y 0 no muestra en pantalla

% nargin es el numero de parametros de la funcion
if nargin<8  % en 6 paramétros display=1
    display=1;
end
NewMAT = Init;% Init=matriz que se completa en cada llamado a función
S=sum2(Zone);
%sum2 suma los elementos de un vector, en este caso de una matriz
% al ejecutar retorna la suma de todos los elementos, siempre da un entero
% posiivo, ya que suma el numero de unos, si por ejemplo la matriz zone es
% de 10 x 12, existen 120 celdas, si retorna sum(b:)=106, quiere decir que
% existen 14 ceros, valores desconocidos(0) por completar y 106 unos.
ii=0;
ErrFrob = inf;
% N - number of iterations=10, tambien puede ser 100
% tol - Frobenius norm tolerance inicialmente es tol = 1e-7;

%con_desv2=con_desv./(S); % en seno realiza 100 operaciones  
%con_desv2=con_desv./(2*S); % en seno realiza 101 operaciones  
%con_desv2=con_desv./(S/2); %en seno realiza 246 operaciones
%con_desv2=con_desv./(S/4); %en seno realiza 306  operaciones, se parece más ala señal
%con_desv2=con_desv./(S/8); %en seno realiza   306 operaciones, se parece más a la señal
con_desv2=con_desv./(S/16); %en seno realiza  294 operaciones, se parece más a la señal

%con_desv2=con_desv./(S/32); %en seno realiza 241 operaciones, se parece más a la señal


% la anterior linea es fundamental ya que se realiza un ajuste ,  esta es
% la adición al algoritmo hecha por olger
%con_desv2=zeros(size(con_desv)); % si no adiciona nada realiza 18 operaciones, pero para nada es la señal que se quiere reconstruir
cal_mat2=0;
while (ii<=N) && (ErrFrob > tol)
    % se ejecuta siempre que las iteraciones sean menores a N(10 o 100)
    % y además siempre que el Errfrob(error frobenius)> tol (frobenius
    % tolerancia)1e-7
    %NewMAT = Init;% Init=matriz que se completa en cada llamado a función
   
    NewMAT = (NewMAT.*(1-Zone))+(Zone.*A)+(con_desv2.*(1-Zone));
    %NewMAT = (NewMAT.*(1-Zone))+(Zone.*A)
    %adiciona olger a lo anterior: +(con_desv2.*(1-Zone));
    % para reducir el numero de operaciones  al encontrar la forma nuclear
    
    
    % 1-Zone,cambia 1 por cero's y viceversa. 0 = perdidos
    % newMat, matriz que se completa en cada iteración
    % A=matriz de origen para completar
    NewMAT = FindNuclearNormApprox(NewMAT, lambda);
    
       
     ErrFrob=sqrt(sum2((Zone.*(NewMAT-A)).^2))/S; %
    %sqrt(sum2((Zone.*(NewMAT-A)).^2) es la norma de frobenius
    % S es el numero de unos en la matriz B=Zone
     fprintf(' Calcula de NewMAT %d \n',ii);
     fprintf('ii=%d, ErrFrob=%g , tolerancia=%g\n',ii,ErrFrob,tol); 
    if display
        fprintf('ii=%d, ErrFrob=%g \n',ii,ErrFrob); 
    end;
    ii=ii+1;
    cal_mat2=ii; % se adiciona para saber el numero de recalculos de la matriz NewMat
    
    Error(ii) = ErrFrob;
    if ii>2
        Diff = Error(ii)-Error(ii-1);
        % cuando el  ErrFrob se parece mucho al anterior, la diferemcia es 0 
        if abs(Diff)<tol/10% tol - Frobenius norm tolerance inicialmente es tol = 1e-7;
            % tol/1000 el valor normal, pero se cambio por 10
            if display
                fprintf('Diff demasiado pequeño cuando el  ErrFrob se parece mucho al anterior \n');
            end
            fprintf('Diff demasiado pequeño cuando el  ErrFrob se parece mucho al anterior \n');
           
            return;
        end
    end
    %fprintf('\n');
end

end



