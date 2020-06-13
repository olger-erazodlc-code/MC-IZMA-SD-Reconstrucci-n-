function B=FindNuclearNormApprox(A, lambda)
% Finds the best LS approximation X, to a givem matrix M, such that
% ||X||<=lambda. The norm is the nuclear norm (sum of singular values).
    [u,s,v]=svd(A,'econ');
    sig = diag(s);
    % obtiene los elementos de la diagonal de una matriz
    % sum suma los  valores por columna de una matriz y el resultado lo
    % genera un vector de 1 fila y n columnas
    % sum aplicado a un vector suma todas las celdas y genera un solo valor
        if sum(sig)<lambda
        fprintf('es verdad sum(sig)<lambda ---> %g < %g \n',sum(sig),lambda);
        % 
        B=A;
        return; % retorna
        end
    n = length(sig);% longitud de los elementos de la diagonal de la matriz
    A = ones(1,n);% crea una matriz de unos, de una fila y  n columnas.
    % A es la matriz de coeficientes lineales de restricción ejemplo:
    % 2x+y<=3, es decir los numeros a cargar son 2   1  por cada restriccion
    
    %A = [ones(1,n-2) 100000 100000];
    b = lambda; % vector constantes en las restricciones lineales
    H = eye(n); % crea la matriz idendidad, 1 en la diagonal y ceros en el resto.
    % Tamaño n x n.por lo general H es la matriz de los coeficientes de terminos
    % cuadráticos
    f = -sig(:);% hace negativo cada elemento 
    % ademas convierte la matriz en una sola columna. f contiene los coeficientes
    % de los terminos lineales.
    lb = zeros(n,1);% crea una matriz de ceros, con n filas y l columnas
    % limites inferiores de las variables, en este caso 0 
    x = quadprog(H,f,A,b,[],[],lb,[],[],optimset('Display','off'));
    %x = qpas(H,f,[],[],A,b,lb,[],0);
    B = u*diag(x)*v';
    % lo  anterior multiplica el vector propio unitario singular por la izquierrda
    % * diagonal de los valores obtenidos en  x que corrresponden a los valores singulares
    % * el vector propio unitario singular por la derecha trasnpuesto
end

    
    
