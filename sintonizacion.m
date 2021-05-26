%%%% SINTONIZACIÓN DE PARAMETROS Y OBTENCIÓN DE GR�?FICAS%%%%%
clear;
clc;
clf;
warning('off');

s=tf('s');

%Parametros del motor
K=2652.28;
p=64.986;
G=K/(s*(s+p));


%% Curvas de Mp
%zeta de 0.2 a 1
v=0.02;

Mp_des=1.115;
beta2=0.5;
betas=0:0.01:40;
zetas=0.1:0.1:1;
L=length(betas);

figure(1);
Legend=cell(length(zetas),1);;
figure(1);
hold on;
n=0;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_param(x,t,v);
        tmp=[tmp Mp];
        
        
    end
    
    n=n+1;    
    
    [beta_inf,beta_sup]=get_betas(tmp,betas,Mp_des);
    Legend{n}=strcat('\zeta=',num2str(zetas(n)),', \beta=(',num2str(beta_inf),',',num2str(beta_sup),')');
    
    plot(betas,tmp,'LineWidth',1.5);
    
    
    
end

title('Curvas de Mp para varios \zeta')
plot(betas,ones(L,1)*Mp_des,'k--','LineWidth',1);
Legend{length(zetas)+1}=strcat('M{p}');
legend(Legend,'FontSize',10);
xlabel('\beta');
ylabel('M{p}');
hold off;
saveas(gcf,'img/betas1.png')
saveas(gcf,'img/betas1.fig')

%% Curvas de Mp
%zeta de 1,2 a 2
v=0.02;

Mp_des=1.115;
beta2=0.5;
betas=0:0.02:40;
zetas=1.1:0.1:2;
L=length(betas);
Legend=cell(length(zetas),1);
figure(2);
hold on;
n=0;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_param(x,t,v);
        tmp=[tmp Mp];
        
        
    end
    
    n=n+1;    
    
    [beta_inf,beta_sup]=get_betas(tmp,betas,Mp_des);
    Legend{n}=strcat('\zeta=',num2str(zetas(n)),', \beta=(',num2str(beta_inf),',',num2str(beta_sup),')');
    
    plot(betas,tmp,'LineWidth',1.5);
    
    
    
end

title('Curvas de Mp para varios \zeta')
plot(betas,ones(L,1)*Mp_des,'k--','LineWidth',1);
Legend{length(zetas)+1}=strcat('M_{p}Max');
legend(Legend,'FontSize',10);
xlabel('\beta');
ylabel('M{p}');
hold off;
saveas(gcf,'img/betas2.png')
saveas(gcf,'img/betas2.fig')

%% Obtencion de betas

v=0.02;

Mp_des=1.115;
beta2=0.5;
betas=0:0.02:36;
zetas=0.1:0.1:2;
L=length(betas);

beta_xls=[];
n=0;

Mps=[];
figure(1);
hold on;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_param(x,t,v);
        tmp=[tmp Mp];
        
    end
    
    n=n+1;
   [beta_inf,beta_sup]=get_betas(tmp,betas,Mp_des);
   beta_xls(n,1)=zeta;
   beta_xls(n,2)=beta_inf;
   beta_xls(n,3)=beta_sup;
   
   
    
end
filename='betas.xls';
xlswrite(filename,beta_xls);

%% Calculo de ts
%para zetas entre 0.2 y 1

v=0.02;



Legend=cell(6,1);
betas2=0.1:0.02:30;
zetas=0.2:0.2:1;
L=length(betas2);
ts_max=0.5;

n=0;
figure(3)

index=0;
index2=0;
hold on;


tuplas_xls=[];

for i=1:(length(zetas))
    
    n=n+1;
   
        tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,2);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;
            for beta2=betas2
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_param(x,t,v);
                    tmp=[tmp ts];
              
            end
            
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;

           
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
        
            tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,3);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;         
            for beta2=betas2
                if zeta==1
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_param(x,t,v);
                    tmp=[tmp ts];
                else
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_param(x,t,v);
                    tmp=[tmp ts];
                end
            

            end
            
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;
            
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
    
    
    
    
     
end
 
title('Curvas de t{s}')
plot(betas2,ones(L,1)*ts_max,'k--','LineWidth',1.2);
Legend{index+1}=strcat('t{s}');
legend(Legend,'FontSize',10);
xlabel('\beta{2}');
ylabel('t{s}');
hold off;
saveas(gcf,'img/curvasTs.fig')
saveas(gcf,'img/curvasTs1.png')

%% Calculo de ts
%para zetas entre 1.2 y 2
%debe ejecutarse justo despues del anterior para que se guarden las tuplas
v=0.02;



Legend=cell(11,1);
betas2=0.1:0.1:40;
zetas=1.2:0.2:2;
L=length(betas2);
ts_max=0.5;

n=5;
figure(4)
Mps=[];
index=0;
hold on;
index2=5;
for i=1:(length(zetas))
    
    n=n+1;
   
        tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,2);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;
            for beta2=betas2
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_param(x,t,v);
                    tmp=[tmp ts];
                    
              
            end
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;

           
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
            tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,3);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;
          
            for beta2=betas2
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_param(x,t,v);
                    tmp=[tmp ts];
             

            end
            
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;
            
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
        
     
end
        filename2='tuplas.xls';         %zetas/betas/betas2
        xlswrite(filename2,tupla_xls);  %0.2  /8.3  /18.6
                                        
 
title('Curvas de t{s}')
plot(betas2,ones(L,1)*ts_max,'k--','LineWidth',1.2);
Legend{index+1}=strcat('t_{s}max');
legend(Legend,'FontSize',18);
xlabel('\beta{2}');
ylabel('t{s}');
hold off;
saveas(gcf,'img/curvasTs2.fig')
saveas(gcf,'img/curvasTs2.png')


%%
%respuesta al escalon
zeta=1.2;
beta=7.8;
beta2=21.4;
figure(6)
    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_param(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)));
    [x,t]=step(H);
    plot(t,x,'LineWidth',1.5);


