clc;clear
K=2652.28;
reductora=23;
Km=K/reductora;
p=64.986;
load tuplas.csv
beta=tuplas(3,2);
beta2=tuplas(3,3)/4;
zeta=tuplas(3,1);
T=5e-3;


[tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,Km,beta,beta2,zeta);



Ki=Kp*T/tau_i;
Kd1=Kp*tau_d1/T;

Kd2=Kp*tau_d2/T;

%%imprimir para copiar y pegar
disp(num2str(Kp,4))
disp(num2str(Ki,10))
disp(num2str(Kd1,10))
disp(num2str(Kd2,10))
%%
disp(num2str(Kp,4))
disp(num2str(tau_d1,4))
disp(num2str(tau_d2,4))
disp(num2str(tau_i,4))
%%
for n=1:15
    beta=tuplas(n,2);
    beta2=tuplas(n,3);
    zeta=tuplas(n,1);
    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,Km,beta,beta2,zeta);
    epsilon = p+K*Kp*tau_d1;
    
    if Kp*tau_d1 <= -p/K
        a=0;
    else
        a=1;
    end
    if epsilon*tau_i>1
        b=1;
    else
        b=0;
    end

    if Kp < 0
        c=0;
    else
        c=1;
    end

    if(c==1&&b==1&&a==1)
        disp(strcat('zeta=',num2str(zeta),', beta=',num2str(beta),' beta_{2}=',num2str(beta2)));
        disp(strcat('estable a= ',num2str(a),' b= ',num2str(b),' c= ',num2str(c)));
    else
        strcat('zeta=',num2str(zeta),', beta=',num2str(beta),' beta_{2}=',num2str(beta2));
        disp(strcat('no estable a= ',num2str(a),' b= ',num2str(b),' c= ',num2str(c)));
        num2
    end
end