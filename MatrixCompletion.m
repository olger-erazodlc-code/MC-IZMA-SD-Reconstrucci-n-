function [NewMAT, ier, cal_mat] = MatrixCompletion(A, Zone,N, Mode, lambda_tol, tol, con_desv, display, Weights )
% Inputs:
% A - initial matrix and matrix to complete as determined by Zone
% N - number of iterations
% Mode - Norm to minimize. Can be "spectral", "nuclear" or
% "NuclearWeighted" (the latter is not really a norm, so global convergence
% is not guaranteed
% lambda_tol - Tolerance for the minimal norm found (smaller than
% theoretical minimum+lambda_tol) inicialmente es lamnbda_tol = 10;
% tol - Frobenius norm tolerance incialmente es tol = 1e-7;
% Weights(pesos) - Give different weights, should be ascending and relevant only
% to "WeightedNuclear" mode.
%
% Output:
% NewMAT - Matrix after completion
% ier - Success/Fail indicator - 0 - No error. =Exito/fracaso
%
% References: 
% G. Shabat, A. Averbuch "Interest Zone Matrix Approximation", Electronic
% Jounal of Linear Algebra, Vol. 23(1), 2012
%
% nargin es el numero de parametros de la funcion}
%fprintf('\nparametros de la funcion:');
%nargin
if nargin<9
    % es 8 si es con norma ponderada
    % es 7 si la norma es nuclear y la norma es espectral, display
    %puede ser 0 ó 1.
    Weights=ones(1,min(size(A)));
    % para la matriz A con entradas perdidas puede ser 3, ya 
    % que tiene 3 filas y 10 columnas
elseif nargin<8
    % es 6 si la norma es nuclear y la norma es espectral, 
    % se asume con display=1
    display=1;
    Weights=ones(1,min(size(A)));
    % crea una matriz de 1's de una fila y el numero de columnas es el mínimo de
    % (#filas y #columnas de A).
end

ier = 0;
min_lambda = 0;
if strcmp(Mode,'nuclear')
    max_lambda = sum(svd(A))*1.1; % A es a matriz con datos perdidos
    % svd=singular valor descomposition
    % svd(A) returns the singular values of matrix A in descending order.
    % es decir que corresponde a S=matriz diagonal de valores singulares 
    % sum realiza la sumatoria de los valores singulares listados en orden
    % descendenete y después se multiplica por un coeficiente-constante C
elseif strcmp(Mode,'spectral')
    max_lambda = max(svd(A))*1.2; % Factor to ease convergence.
elseif strcmp(Mode,'NuclearWeighted')
    s=svd(A);
    max_lambda = Weights*s;
else
    fprintf('Undefined mode. quitting. \n');
    return;
end
%tol = 0.00001;
NewMAT = A; % asigna la matriz incompleta a completar
lambda = inf; % inf retorna un numero positivo infinito
lambda_prev = 0;
%lambda_tol = 1;
err = inf; % incialmente el err le asigna infinito
Counter = 0;
Converge=0;
%NewMAT=A.*Zone;
% cuando inicia lamnbda_tol = 10; tol = 1e-7; err=infinito;lambda_prev=0;lamnbda_tol = 10;
% while ((infinito > 1e-7 ) || (abs(infinito-0) > 10))
fprintf('Matriz que tiene los valores semilla con desviaion estandar \n');
con_desv
cal_mat=0;% se adiciona para saber el numero de operaciones acumuladas para el recalculo de la matriz, al aproximar
S=sum2(Zone);
con_desv2=con_desv./(S*15);
while ((err > tol) || (abs(lambda-lambda_prev) > lambda_tol))
    Counter=Counter+1;
    lambda_prev = lambda; % la primera vez es infinito, la segunda vez 22
    lambda=(min_lambda+max_lambda)/2;% la primera vez min_lambda=0 y max_lambda=45 actualiza a 22 
    if strcmp(Mode,'nuclear')
        NewMAT =NewMAT+(con_desv2.*(1-Zone));
        [NewMAT, Error,cal_mat2] = MatApproxNuclear(A.*Zone, NewMAT, Zone, lambda, N, tol, con_desv, display); % Was A.*Zone instead NewMAT in the input
        % lo anterior es una modificación
        %A.*Zone es la matriz imcompleta original, multiplicada por la matriz de 0 y
        %1, donde 0 representa la perdidad de datos
        % NewMat es la matriz a completar, la cual en cada ejecución se
        % actualiza con el resultado MatApproxNuclear, es decir que se va
        % completando en cada ejecución.
        % Zone=matriz de ceros y unos, ceros son datos perdidos =antes en completion_demo se llama B
        % lambda 
        % N es el número de iteraciones 
        % tolerancia norma frobenius tol= 1e-7
        % display=1,  si la norma es nuclear o la norma es espectral
        % muestra resultados
        %NewMAT=NewMAT+ ((1-Zone).*con_desv)
        cal_mat=cal_mat+cal_mat2;
    elseif strcmp(Mode,'spectral');
        [NewMAT, Error] = MatApproxSpectral(A.*Zone, NewMAT, Zone, lambda, N, tol);
    elseif strcmp(Mode, 'NuclearWeighted')
        [NewMAT, Error] = MatApproxNuclearWeighted(A.*Zone, NewMAT, Zone, Weights, lambda, N, tol, display);
    else
        fprintf('Modo indefinido. Abandonando. \n');
        return;
    end
    err = Error(end);% extrae el ultimo valor de la matriz Error
    % es decir el cruce del maximo numero de la fila y de la columna
    % la cual se obtiene por la norma:  ||pX - pM||frobenius
    % pX es la nueva matriz que se completa en cada iteracion= NewMAT
    % pM es la matriz a completar A.*Zone
    if Error(end) > tol  % tol - Frobenius norm tolerance incialmente es tol = 1e-7;
        min_lambda=lambda; % si el error es > a tol, 
    else
        max_lambda=lambda;
    end
    if abs(lambda-lambda_prev) < lambda_tol
        Converge=Converge+1;
    end
    fprintf('_____________________converge %u___________________\n',Converge);
    NewMAT
    fprintf('Hasta el momento, numero de convergencias es: %d \n', Converge);
    if Converge>4
        if display
            fprintf('Parece que el algoritmo no pudo converger. Hay dos opciones: \n');
            fprintf('La norma de la matriz inicial es demasiado pequeña o se necesitan más iteraciones (N mayor) \n');
        end
        ier=1;
        return;
    end
    if display
        fprintf('Error=%g, lambda=%g, lambda-lambda_prev=%g, Counter=%d \n',Error(end),lambda, lambda-lambda_prev,Counter);
    end
     fprintf('Error=%g, lambda=%g, lambda-lambda_prev=%g, Counter=%d \n',Error(end),lambda, lambda-lambda_prev,Counter);
end % final del while

end
