%% --------------------------------------------SENO-------------------------------------------------------------
close all;% cierra figuras
clear all;%Limpia el espacio de trabajo
clc; % limpia 
isimula=0;
var_ruido=0.0; %antes 0.2 para simular y para documento tiende a cero

%debe ir desde 0.1 hasta 4.0
disp('¡Hola olger inicia el proceso de simuacion continua!') %% señal seno a transmitir
while(true)
    cad=  matrix_completion_izma_sd_auto(0,1,var_ruido);
    % tipo_modulacion,senal_mas_ruido,sigma_v
     %   tipo_modulacion=0 seno
    %   tipo_modulacion=1 coseno
    %   tipo_modulacion=2 solo ruido
    %   tipo_modulacion=3 bpsk
    %   tipo_modulacion=4 qpsk
    
    %   senal_mas_ruido=0  solo ruido
    %   senal_mas_ruido=1  señal más ruido
    
    % desviacion estandar 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
         
    addpath('D:\maestria\semestre13\script tesis completion_con_ruido_mod_simula_auto','-end');
    formatOut = 'dd_mm_yyyy';
    fhoraa=datestr(now,formatOut);
    nombre_archivo=strcat('D:\maestria\semestre13\script tesis completion_con_ruido_mod_simula_auto\resultados\seno_resul_simula_',fhoraa,'_0_01.txt');
    fid=fopen(nombre_archivo,'a+'); 
    fprintf(fid,cad); 
    fclose(fid);
    
    % lo siguiente indica los resultados de la simulacion
    isimula=isimula+1;
    str_simula = sprintf('Termino la simulacion número: %d', isimula);
    disp(str_simula);
    disp('___________________________________________________________');
    var_ruido=var_ruido+0.0035; % al termino de las dos mil veces se hace 1
       
    pause(1);
    %cada fila tiene 65 caracteres, vamos a trabajar con 100 registros por 
    % simulacion, por ello la longitud debe ser inferior a 7000 cxaracteres
    fid2=fopen(nombre_archivo);
    contenido=fread(fid);
    fclose(fid2);
    [fila col]=size(contenido);

    if fila>120000
        break;
    end
 break;
end