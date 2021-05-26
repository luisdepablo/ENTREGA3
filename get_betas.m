function [beta_inf,beta_sup]=get_betas(tmp,betas,Mp_des)
beta_inf=0;
beta_sup=0;
L=length(betas);

for i=1:(L-1)
   
   
        if (tmp(i)==Mp_des && tmp(i+1)>Mp_des)||(tmp(i)<=Mp_des && tmp(i+1)>=Mp_des)||(tmp(i)>=Mp_des && tmp(i-1)<=Mp_des)
            beta_inf=betas(i);
            

        elseif(tmp(i)==Mp_des && tmp(i+1)<Mp_des)||(tmp(i)<=Mp_des && tmp(i-1)>=Mp_des)||(tmp(i)>=Mp_des && tmp(i+1)<Mp_des)
            beta_sup=betas(i);
            
        end
   
    
end

