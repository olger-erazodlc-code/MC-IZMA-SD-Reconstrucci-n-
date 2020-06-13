function cadena_res_detec= matrix_completion_izma_sd_auto(tipo_modulacion,senal_mas_ruido,sigma_v)
% Problem Size: n1, n2 - matrix dimensions, r - rank
% el rango indica el numero de filas o columnas independienestes
amplitud=1;
T=linspace(0,1,24)';
seno1=amplitud*sin(2*pi*9*T);
if (tipo_modulacion==0)
    seno1=amplitud*sin(2*pi*9*T);
    
end
if (tipo_modulacion==1)
    seno1=amplitud*cos(2*pi*9*T);
end
if (tipo_modulacion==2)
    ruido=0+sigma_v.*randn(size(seno1));
    seno1=ruido; %xR es la Señal recibida mas el ruido
end
if (tipo_modulacion==3)
    %% BPSK
    display('Binario PSK');
	freq=3;
	T=linspace(0,1,10)';
    for ii = 1:1:10 % EN LA PARTE SUPERIOR INDICA 10 LA LONGITUD
        car1(ii) = sin((2*pi*freq*T(ii))); %CARRIER TO BE TRANSMITTED
        % Portadora a transmitir
       
    end
    bpsk1=car1;
    bpsk0=[-1*car1];
    bpsk=[bpsk1 bpsk0(2:10) bpsk1(2:6) ];
    
    %T21=linspace(0,1,24)';
	n11=linspace(0,1,24)';
	t11=linspace(0,1,1/0.001)'; % es lo mismo que un vector desde 0 hasta 1, que contenga 1000 elementos 
	Tm11=1/24;
	parte1=(1/Tm11)*t11(:,ones(size(n11)));
	parte2=(1/Tm11)*n11(:,ones(size(t11)))';
	parte3=parte1-parte2;
	parte4=sinc(parte3);
	parte55=parte4*bpsk';
	  	
    seno1=bpsk';
end
fprintf('Señal completa de 24 muestras original: \n');
seno1
if (tipo_modulacion==4)
    %% QPSK
    display('Cuadratura PSK');
	freq=3;
	T=linspace(0,1,10)';
    for ii = 1:1:10 % EN LA PARTE SUPERIOR INDICA 10 LA LONGITUD
        car1(ii) = sin((2*pi*freq*T(ii))-135); %CARRIER TO BE TRANSMITTED 00
        car2(ii) = sin((2*pi*freq*T(ii))-45); %CARRIER TO BE TRANSMITTED  01
     end
    qpsk1=car1;
    qpsk0=car2;
    qpsk=[qpsk1 qpsk0(2:10) qpsk1(2:6) ];
    
    %T21=linspace(0,1,24)';
	n11=linspace(0,1,24)';
	t11=linspace(0,1,1/0.001)'; % es lo mismo que un vector desde 0 hasta 1, que contenga 1000 elementos 
	Tm11=1/24;
	parte1=(1/Tm11)*t11(:,ones(size(n11)));
	parte2=(1/Tm11)*n11(:,ones(size(t11)))';
	parte3=parte1-parte2;
	parte4=sinc(parte3);
	parte55=parte4*qpsk';
	  	
    seno1=qpsk';
end

seno_puro=seno1;
if (senal_mas_ruido==0)
% solo ruido
        ruido=0+sigma_v.*randn(size(seno1));
        seno1=ruido; %xR es la Señal recibida mas el ruido
end
if (senal_mas_ruido==1)
        % señal mas ruido
        %sigma_v=2.5;
        ruido=0+sigma_v.*randn(size(seno1));
         % ruido  es una matriz aleatoria del tamaño de xR,
        % con media en 0 y desviacion estandar=sigma_v
        seno1=seno1+ruido; %xR es la Señal recibida mas el ruido
end

%% calculo de SNR

potencia_senal=var(seno_puro)+mean(seno_puro)^2;
potencia_ruido=var(ruido)+mean(ruido)^2;
if (potencia_ruido<potencia_senal)
db_snr=potencia_senal/potencia_ruido;
else
db_snr=-(potencia_ruido/potencia_senal);
end
 str_simula = sprintf('SNR de la señal: %02f', db_snr);
 disp(str_simula);

%seno1=awgn(seno1,30);
% si solo se quiere ruido, entonces:
%seno1=ruido; %xR es la Señal recibida pero solo ruido
% en la operacion anterior se multiplica la matriz aleatoria por la
% desviacion estandar y se le suma a la señal recibida para agregar ruido
%% graficar la señal con ruido
n1=linspace(0,1,24)';
t1=linspace(0,1,1/0.001)'; % es lo mismo que un vector desde 0 hasta 1, que contenga 1000 elementos 
%t=linspace(0,1,1000)
Tm1=1/24;
nn1=linspace(0,1,250);
if (tipo_modulacion==0)
    seno21=amplitud*sin(2*pi*9*nn1);
end
if (tipo_modulacion==1)
    seno21=amplitud*cos(2*pi*9*nn1);
end


%ya1=sinc((1/Tm1)*t1(:,ones(size(n1)))-(1/Tm1)*n1(:,ones(size(t1)))')*seno1';
figure(1);
if (tipo_modulacion==0)
    parte1=(1/Tm1)*t1(:,ones(size(n1)));
    parte2=(1/Tm1)*n1(:,ones(size(t1)))';
    parte3=parte1-parte2;
    parte4=sinc(parte3);
    parte5=parte4*seno1;
    plot(n1,seno1,'or',t1,parte5,'b',nn1,seno21,'.-g');
    grid;
    legend('Muestras con ruido','Señal seno con ruido','Señal sin ruido');
end
if (tipo_modulacion==1)
     parte1=(1/Tm1)*t1(:,ones(size(n1)));
    parte2=(1/Tm1)*n1(:,ones(size(t1)))';
    parte3=parte1-parte2;
    parte4=sinc(parte3);
    parte5=parte4*seno1;
    plot(n1,seno1,'or',t1,parte5,'b',nn1,seno21,'.-g');
    grid;
    legend('Muestras con ruido','Señal coseno con ruido','Señal seno sin ruido');
end
if (tipo_modulacion==2)
    parte1=(1/Tm1)*t1(:,ones(size(n1)));
    parte2=(1/Tm1)*n1(:,ones(size(t1)))';
    parte3=parte1-parte2;
    parte4=sinc(parte3);
    parte5=parte4*seno1;
    plot(n1,seno1,'or',t1,parte5,'b');
    grid;
    legend('Muestras del ruido','Señal del ruido');
end
if (tipo_modulacion==3)
    plot(n11,bpsk,'or',t11,parte55,'b');
    grid;
    legend('Muestras','portadora bpsk 24 muestras');
    title('portadora bpsk a transmitir')
    xlabel('Muestras')
    ylabel('Amplitud')
    grid;
 end
if (tipo_modulacion==4)
    plot(n11,qpsk,'or',t11,parte55,'b');
    grid;
    legend('Muestras','portadora qpsk 24 muestras');
    title('portadora qpsk a transmitir')
    xlabel('Muestras')
    ylabel('Amplitud')
    grid;
 end

xlabel('Tiempo');
ylabel('Amplitud');
axis([0 1 -1.9 1.9]);


%% se continua con la señal

figure(2);
%plot(seno1);%antes para graficar como lineas
stem(seno1);%para graficar puntos
if (tipo_modulacion==0)
title('señal seno1 a transmitir sin perdidas')    
end
if (tipo_modulacion==1)
title('señal coseno1 a transmitir sin perdidas')    
end
if (tipo_modulacion==2)
title('señal ruido a transmitir sin perdidas')    
end
if (tipo_modulacion==3)
title('señal bpsk a transmitir sin perdidas')    
end
if (tipo_modulacion==4)
title('señal qpsk a transmitir sin perdidas')    
end
xlabel('Muestra')
ylabel('Amplitud')

seno22=[seno1(2:24); seno1(2)];
ACOMPLETA=[seno1'; seno22' ;seno1'; seno22';seno1'; seno22' ;seno1'; seno22'];
%% señal que se transmite en este caso una función seno 9 hz

% B = rand(size(A))<0.7; % remove 30% of the entries=0.7
B = [
    1 0 1 0 1 0 1 0 1 0 1 0       1 0 1 0 1 0 1 0 1 0 1 0; 
    0 1 0 1 0 1 0 1 0 1 0 1       0 1 0 1 0 1 0 1 0 1 0 1; 
    1 0 1 0 1 0 1 0 1 0 1 0       1 0 1 0 1 0 1 0 1 0 1 0; 
    0 1 0 1 0 1 0 1 0 1 0 1       0 1 0 1 0 1 0 1 0 1 0 1;
    1 0 1 0 1 0 1 0 1 0 1 0       1 0 1 0 1 0 1 0 1 0 1 0; 
    0 1 0 1 0 1 0 1 0 1 0 1       0 1 0 1 0 1 0 1 0 1 0 1; 
    1 0 1 0 1 0 1 0 1 0 1 0       1 0 1 0 1 0 1 0 1 0 1 0; 
    0 1 0 1 0 1 0 1 0 1 0 1       0 1 0 1 0 1 0 1 0 1 0 1   
    ]; % 50% de perdidas
numero_de_unos=sum2(B);
[filas,columnas] = size(B);
perdidos=(filas*columnas)-numero_de_unos;
fprintf('\nDatos perdidos %g de %g \n',perdidos,filas*columnas);

figure(3);
%plot( ACOMPLETA(1,:).*B(1,:),'red'); %grafica con rectas
stem( ACOMPLETA(1,:).*B(1,:),'red'); %grafica con puntos
if (tipo_modulacion==0)
title('señal seno1 a transmitir con perdidas')    
end
if (tipo_modulacion==1)
title('señal coseno1 a transmitir con perdidas')    
end
if (tipo_modulacion==2)
title('señal ruido a transmitir con perdidas')    
end
if (tipo_modulacion==3)
title('señal BPSK a transmitir con perdidas')    
end
if (tipo_modulacion==4)
title('señal QPSK a transmitir con perdidas')    
end
xlabel('Muestra')
ylabel('Amplitud')

% ahora se procede a sembrar una  semilla ya que es un sinusoidal, por ello
% en los espacios de 0 se cambia por un valor que corrresponde a
% desviación estandar de cada una de las filas, se asigna positivo el valor
% de la desviación si el anterior y el siguiente son positivos, al igual
% que si es positivo uno de ellos  y el otro es cero.
% si el anterior y el siguiente son negativos se procede a cargar la
% semilla con la desviación estandar pero negativa
a_integrar_semilla=(ACOMPLETA(1,:).*B(1,:));
[semilla1]=integre_semilla_desv_stan4(a_integrar_semilla,1); 
fprintf('antes de cargar semilla: \n');
a_integrar_semilla
fprintf('Después de cargar semilla: \n');
semilla1

% ultimo parametro indica con 1 si inicia con 1 
a_integrar_semilla=(ACOMPLETA(2,:).*B(2,:));
[semilla2]=integre_semilla_desv_stan4(a_integrar_semilla,0); 


figure(4);
%plot(semilla1); antes grafica con lineas
stem(semilla1); %grafica con puntos
if (tipo_modulacion==0)
title('señal seno con datos perdidos y desv estandar')
end
if (tipo_modulacion==1)
title('señal coseno con datos perdidos y desv estandar')
end
if (tipo_modulacion==2)
title('Ruido con datos perdidos y desv estandar')
end
if (tipo_modulacion==3)
title('BPSK con datos perdidos y desv estandar')
end
if (tipo_modulacion==4)
title('QPSK con datos perdidos y desv estandar')
end
xlabel('Muestra')
ylabel('Amplitud')


ASD=[semilla1; semilla2;semilla1;semilla2;semilla1; semilla2;semilla1;semilla2];
A=ACOMPLETA;
%A=ACOMPLETA.*B;
rango=rank(A);
fprintf('\nMatriz tiene rango: %g \n',rango);
%B es una matriz con unos y ceros aleatorios, 0's representa datos perdidos
lamnbda_tol = 10;%Valor de tolerancia de la norma nuclear / valor espectral - valor mínimo.
tol = 1e-7;%Tolerancia en las entradas conocidas.
N = 50;% originalmente es 100 , esto lo adiciona olger para reducir tiempo

fprintf('Completar usando la norma nuclear de minimización... \n');
[CompletedMat, ier,calculos_matrices] = MatrixCompletion(A.*B, B,N, 'nuclear', lamnbda_tol, tol, ASD,0);

% envia como parametro la matriz con datos perdidos =A.*B
% matriz de ceros y unos, ceros son datos perdidos =B
% numero de iteraciones =N, maximo 100
% norma para minimizar 'nuclear'
% tolerancia para la minima norma encontrada=lamnbda_tol=10
% tolerancia norma frobenius tol= 1e-7
% A es la matriz que tiene desviación estandar donde no existian muestras.

nor_nuc_ini=sum(svd(A.*B));
fprintf('\nMatriz corrupta norma nuclear (inicial): %g \n',nor_nuc_ini);

nor_nuc_fin=sum(svd(CompletedMat));
fprintf('Matriz restaurada norma nuclear (final): %g \n',nor_nuc_fin);

err_med_cua=sqrt(sum2((CompletedMat-A).*B)/sum(B(:)));
fprintf('MSE - error medio cuadrado en entradas conocidas: %g \n',err_med_cua);
fprintf('\nMatriz que se completa:');
CompletedMat

figure(5);
%plot(CompletedMat(1,:)); grafica con lineas
stem(CompletedMat(1,:));
if (tipo_modulacion==0)
title('señal seno recuperada con Norma Nuclear')
end
if (tipo_modulacion==1)
title('señal coseno recuperada con Norma Nuclear')
end
if (tipo_modulacion==2)
title('señal ruido recuperada con Norma Nuclear')
end
if (tipo_modulacion==3)
title('señal bpsk recuperada con Norma Nuclear')
end
if (tipo_modulacion==4)
title('señal qpsk recuperada con Norma Nuclear')
end
xlabel('muestra')
ylabel('Amplitud')


% promedio_senal=CompletedMat(1,:)+CompletedMat(2,:)+CompletedMat(3,:)+CompletedMat(4,:)+CompletedMat(5,:)+CompletedMat(6,:)+CompletedMat(7,:)+CompletedMat(8,:);

promedio_senal=CompletedMat(1,:)+CompletedMat(3,:);%+CompletedMat(3,:)+CompletedMat(4,:);
promedio_senal=promedio_senal./2% se divide entre el exponente del 2, que es igual al numero de columnas.
% ejemplo: 8 columnas, luego 2^3 =8, divide entre 3
figure(6);
plot(promedio_senal);
if (tipo_modulacion==0)
title('señal seno recuperada como promedio')
end
if (tipo_modulacion==1)
title('señal coseno recuperada como promedio')
end
if (tipo_modulacion==2)
title('señal ruido recuperada como promedio')
end
if (tipo_modulacion==3)
title('señal bpsk recuperada como promedio')
end
if (tipo_modulacion==4)
title('señal qpsk recuperada como promedio')
end
xlabel('Tiempo')
ylabel('Amplitud')


n=linspace(0,1,24)';
t=linspace(0,1,1/0.001)'; % es lo mismo que un vector desde 0 hasta 1, que contenga 1000 elementos 
%t=linspace(0,1,1000)
Tm=1/24;
nn=linspace(0,1,250);
if (tipo_modulacion==0)
seno2=amplitud*sin(2*pi*9*nn);
end
if (tipo_modulacion==1)
seno2=amplitud*cos(2*pi*9*nn);
end


ya=sinc((1/Tm)*t(:,ones(size(n)))-(1/Tm)*n(:,ones(size(t)))')*promedio_senal';


figure(7);
if (tipo_modulacion==0)
    plot(n,seno1,'or',n,semilla1,'*c',t,ya,'.-g',nn,seno2,'b',n,promedio_senal,'*y');
    grid;
    legend('Muestras','Semilla SD','Seno reconstruida','Seno Original','IZMA-SD');
 end
if (tipo_modulacion==1)
    plot(n,seno1,'or',n,semilla1,'*c',t,ya,'.-g',nn,seno2,'b',n,promedio_senal,'*y');
    grid;
    legend('Muestras','Semilla SD','Coseno reconstruida','Coseno Original','IZMA-SD');
 end
if (tipo_modulacion==2)
    plot(n,seno1,'or',n,semilla1,'*c',t,ya,'.-g',t1,parte5,'b',n,promedio_senal,'*y');
    grid;
    legend('Muestras','Semilla SD','Ruido reconstruida','Ruido Original','IZMA-SD');
end
 if (tipo_modulacion==3)
    plot(n,seno1,'or',n,semilla1,'*c',t,ya,'.-g',t11,parte55,'b',n,promedio_senal,'*y');
    grid;
    legend('Muestras','Semilla SD','BPSK reconstruida','BPSK Original','IZMA_SD');
 end
 if (tipo_modulacion==4)
    plot(n,seno1,'or',n,semilla1,'*c',t,ya,'.-g',t11,parte55,'b',n,promedio_senal,'*y');
    grid;
    legend('Muestras','Semilla SD','QPSK reconstruida','QPSK Original','IZMA-SD');
 end
xlabel('Tiempo');
ylabel('Amplitud');
axis([0 1 -1.9 1.9]);

%% ahora se calcula que tan semejantes son las selñales

%ya= tiene  1000 elementos, ahora debe tener 250 muestras, ya que esta es
   %la señal reconstruida para así compararla con la original
   %seno2= tiene 250 muestras de la señal original
 [fila,colu] = size(ya) ;%devuelve el número de filas y columnas cuando A es una matriz.
 fprintf('\nfilas de ya: %g    y columnas de ya: %g \n',fila,colu);
 %ahora debemos construir una matriz de 250 elementos
 if (tipo_modulacion==0) || (tipo_modulacion==1)
     reconstruida_a_comparar=zeros(1,250);
     txx=0;
    
     for ii=1:fila 
            residuo= mod(ii, 4);
            if residuo==0
                 txx=txx+1;
                recorte=ya(ii,1);
                %fprintf('%d\n',recorte);
               reconstruida_a_comparar(1,txx)=recorte';
            end 
     end
     max_senal=max(reconstruida_a_comparar(1,1:250));
     escala_reconstruida_a_comparar=reconstruida_a_comparar(1,1:250)./(max_senal); 

     max_senal=max(seno2(1,1:250));
     escala_seno2=seno2(1,1:250)./(max_senal); 
     [res0,res1]= xcorr( escala_reconstruida_a_comparar,escala_seno2);
 end     
  
 if (tipo_modulacion==2)
      parte5x=parte5';
      max_senal=max(parte5x(1,1:1000));
      escala_seno2=parte5x(1,1:1000)./(max_senal); 
      max_senal=max(ya(1:1000,1));
      escala_ya= ya(1:1000,1)./(max_senal); 
      
          
     reconstruida_a_comparar_250=zeros(1,250);
     escala_seno_250=zeros(1,250);
     txx=0;
      for ii=1:fila 
            residuo= mod(ii, 4);
            if residuo==0
                 txx=txx+1;
                recorte=escala_ya(ii,1);
                %fprintf('%d\n',recorte);
              reconstruida_a_comparar_250(1,txx)=recorte';
              recorte= escala_seno2(1,ii);
                %fprintf('%d\n',recorte);
               escala_seno_250(1,txx)=recorte';
            end 
     end
      
     [res0,res1]= xcorr(reconstruida_a_comparar_250,escala_seno_250);
    %  [res0,res1]= xcorr(escala_seno_250,escala_seno_250);
    %[res0,res1]= xcorr(reconstruida_a_comparar_250,reconstruida_a_comparar_250);
 end
 
 if (tipo_modulacion==3) || (tipo_modulacion==4)
        %[fila2w,colu2w] = size(parte55) ;
        %fprintf('\nfilas de parte55: %d    y columnas de ya: %d \n',fila2w,colu2w);
       parte555=parte55';
       max_senal=max(parte555(1,1:1000));
       escala_seno2=parte555(1,1:1000)./(max_senal); 
       
        max_senal=max(ya(1:1000,1));
       escala_ya= ya(1:1000,1)./(max_senal); 
      
        reconstruida_a_comparar_250=zeros(1,250);
        escala_seno_250=zeros(1,250);
     txx=0;
      for ii=1:fila 
            residuo= mod(ii, 4);
            if residuo==0
                 txx=txx+1;
                recorte=escala_ya(ii,1);
                %fprintf('%d\n',recorte);
                reconstruida_a_comparar_250(1,txx)=recorte';
                recorte= escala_seno2(1,ii);
                %fprintf('%d\n',recorte);
               escala_seno_250(1,txx)=recorte';
            end 
     end
      
     [res0,res1]= xcorr(reconstruida_a_comparar_250,escala_seno_250); 
 end    
 

 r1=sum(res0);
 fprintf('\nsuma de la correlacion cruzada entre la original y la recuperada: %g \n',r1);
 [fila2,colu2] = size(res0) ;%devuelve el número de filas y columnas cuando A es una matriz.
 semejanza=0;
 max_senal2=max(res1(1,1:colu2));
 res11=res1(1,1:colu2)./(max_senal2); 
 
 max_senal2=max(res0(1,1:colu2));
 res00=res0(1,1:colu2)./(max_senal2); 
 
 for ii=1:colu2
     recorte=res00(1,ii).*res11(1,ii);
     semejanza=semejanza+recorte;
 end
 semejanza=abs(semejanza);
 % si semejanza es igual a cero, slas señales se parecen
 fprintf('valor de semejanza %0.4f ---> si esta cerca de 0 se parece \n',semejanza);
    
 figure(8);
 stem(res1,res0);
    

%% ______________ ya termino de simular ________________
   
 
    a=1; 
    %format0:
    % $modulacion:bpsk-qpsk   *senal_mas_ruido     _sigma_vl      
    formatOut = 'dd/mm/yyyy HH:MM:ss';
    fhora=datestr(now,formatOut);
    if (tipo_modulacion==0)
        %str_simula = sprintf('$seno   *%s    _dsr%0.2f   ¿snr%0.2f       ?nni%0.2f   ?nnf%0.2f  #emc%0.2f     (ope_cal_mat%d     +%s\r\n',senal_mas_ruido,sigma_v,db, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,fhora);
        str_simula = sprintf('$seno*%d_dsr%0.3f>snr%0.2f?nni%0.2f?nnf%0.2f#emc%0.12f(ocm%d)sem%0.4f+%s\r\n',senal_mas_ruido,sigma_v,db_snr, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,semejanza,fhora);
    end
    if (tipo_modulacion==1)
        str_simula = sprintf('$cose*%d_dsr%0.3f>snr%0.2f?nni%0.2f?nnf%0.2f#emc%0.12f(ocm%d)sem%0.4f+%s\r\n',senal_mas_ruido,sigma_v,db_snr, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,semejanza,fhora);
    end
    if (tipo_modulacion==2)
        str_simula = sprintf('$ruid*%d_dsr%0.3f>snr%0.2f?nni%0.2f?nnf%0.2f#emc%0.12f(ocm%d)sem%0.4f+%s\r\n',senal_mas_ruido,sigma_v,db_snr, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,semejanza,fhora);
    end
    if (tipo_modulacion==3)
        str_simula = sprintf('$bpsk*%d_dsr%0.3f>snr%0.2f?nni%0.2f?nnf%0.2f#emc%0.12f(ocm%d)sem%0.4f+%s\r\n',senal_mas_ruido,sigma_v,db_snr, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,semejanza,fhora);
    end
    if (tipo_modulacion==4)
        str_simula = sprintf('$qpsk*%d_dsr%0.3f>snr%0.2f?nni%0.2f?nnf%0.2f#emc%0.12f(ocm%d)sem%0.4f+%s\r\n',senal_mas_ruido,sigma_v,db_snr, nor_nuc_ini,nor_nuc_fin,err_med_cua,calculos_matrices,semejanza,fhora);
    end
    display(str_simula);
    cadena_res_detec=str_simula;
    
    
end


