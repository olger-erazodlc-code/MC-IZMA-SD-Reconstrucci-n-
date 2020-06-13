function [semilla2]=integre_semilla_desv_stan4(mat_incon, inicio2);
semilla=mat_incon;
[filas_sem,columnas_sem]=size(mat_incon);
desv1=std2(semilla);
fprintf('total de columnas de mat_incon(integre_semi): %d y el valor de inicio es: %d \n', columnas_sem, inicio2);
for i=1:columnas_sem

    if (inicio2==1)
    %1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0;         
        res = mod(i,2);
        if (res==0) && (i<columnas_sem)
           anterior=semilla(1,i-1);
           proximo=semilla(1,i+1);
            if (anterior<0) && (proximo<0)
                semilla(1,i)=1*desv1;
              
            end
            if (anterior>=0) && (proximo>=0)
                semilla(1,i)=-desv1;
            end
            
            
            if (anterior>=0) && (proximo<0) % 1    y  -2  
                 xy=anterior+proximo;
                 xy=xy/2;
                 semilla(1,i)=-1*xy;
            end
            
             if (anterior<0) && (proximo>=0)
                xy=anterior+proximo;
                 xy=xy/2;
                 semilla(1,i)=-1*xy;
             end
                      
             aux1=round(anterior*10^2)/10^2; 
             aux2=round(proximo*10^2)/10^2; 
             aux3=abs(aux1)-abs(aux2);
             if (abs(aux3)<=0.2)
                if (anterior>=0) && (proximo>=0)
                    semilla(1,i)=-1*(abs(aux3));   
                else
                    semilla(1,i)=abs(aux3);   
                end   
             end
            
        end
        
          if (i==columnas_sem)&& (res==0)
            anterior=semilla(1,i-1);
            if (anterior<0)
                %semilla(1,i)=-1*desv1/2;
                semilla(1,i)=1*desv1/1;
            else
                %semilla(1,i)=desv1/2;
                semilla(1,i)=-1*desv1/1;
            end
          end
          
    end % fin de incio2==1
    
   %% si inicia con 0
    if (inicio2==0)
    %0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1;         
    res = mod(i,2);
    
    if (res==1) && (i<columnas_sem) && (i>1)
       anterior=semilla(1,i-1);
       proximo=semilla(1,i+1);
        if (anterior<0) && (proximo<0)
            semilla(1,i)=1*desv1;
           
        end
        if (anterior>=0) && (proximo>=0)
            semilla(1,i)=-1*desv1;
          
        end
        if (anterior>=0) && (proximo<0)
            xy=anterior+proximo;
            xy=xy/2;
            semilla(1,i)=-1*xy;
        end
         if (anterior<0) && (proximo>=0)
            xy=anterior+proximo;
            xy=xy/2;
            semilla(1,i)=-1*xy;
         end
         
            aux1=round(anterior*10^2)/10^2; 
             aux2=round(proximo*10^2)/10^2; 
             aux3=abs(aux1)-abs(aux2);
             if (abs(aux3)<=0.2)
                if (anterior>=0) && (proximo>=0)
                    semilla(1,i)=-1*(abs(aux3));   
                else
                    semilla(1,i)=abs(aux3);   
                end   
             end
         
    end
            
    
    if (i==1)
        proximo=semilla(1,i+1);
        if (proximo<0)
            %semilla(1,i)=1*desv1/2;
            semilla(1,i)=1*desv1/1;
        else
            %semilla(1,i)=-1*desv1/2;
            semilla(1,i)=-1*desv1/1;
        end
        
    end
    
    if (i==columnas_sem)&&(res==1)
        anterior=semilla(1,i-1);
        if (anterior<0)
            %semilla(1,i)=-1*desv1/2;
            semilla(1,i)=1*desv1/1;
        else
            %semilla(1,i)=desv1/2;
            semilla(1,i)=-desv1/1;
        end
    end
    
    end

end
   semilla2=semilla;  
end






