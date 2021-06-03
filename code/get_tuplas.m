function [tupla]=get_tuplas(tmp,ts_max,betas2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L=length(betas2);

for i=3:(L-2)
   
   
        if (tmp(i)==ts_max && tmp(i+1)>ts_max)||(tmp(i)<=ts_max && tmp(i+1)>=ts_max)||(tmp(i)>=ts_max && tmp(i-1)<=ts_max)
            tupla=betas2(i);
        end

%         elseif(tmp(i)==ts_max && tmp(i+1)<ts_max)||(tmp(i)<=ts_max && tmp(i-1)>=ts_max)||(tmp(i)>=ts_max && tmp(i+1)<ts_max)
%             tupla=betas1(i);
%             
%         end
%    
    
end

end




