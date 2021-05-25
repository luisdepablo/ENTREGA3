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
Mps=[];

betas=0:0.01:40;
zetas=0.2:0.2:2;
L=length(betas);

figure(1);

Legend=cell(length(zetas),1);

Mps=[];
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
legend(Legend,'FontSize',20);
xlabel('\beta');
ylabel('M{p}');
hold off;
saveas(gcf,'img/betas1.png')

%% Curvas de Mp
%zeta de 1,2 a 2
v=0.02;

Mp_des=1.115;
beta2=0.5;
Mps=[];

betas=0:0.01:40;
zetas=1.2:0.2:2;
L=length(betas);

Legend=cell(length(zetas),1);

Mps=[];
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
Legend{length(zetas)+1}=strcat('M{p}');
legend(Legend,'FontSize',20);
xlabel('\beta');
ylabel('M{p}');
hold off;
saveas(gcf,'img/betas2.png')

%% Obtención de betas

v=0.02;

Mp_des=1.115;
beta2=0.5;
Mps=[];
leyenda=[];
betas=0:0.01:40;
zetas=0.2:0.2:2;
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

