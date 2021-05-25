function [Mp,ts,tp,tr]=get_param(x,t,v)
    l=length(x);

    Mp=0;
    tr=t(l);
    tp=0;
    ts=t(l);

    for i=1:l
       if x(i)>Mp
           Mp=x(i);
           tp=t(i);
       end

       if x(i)>1+v
           ts=t(i);
       elseif x(i)<1-v
           ts=t(i);

       end

    end

    for i=l
        if x(i)>1
            tr=t(i);
            break
        end
    end