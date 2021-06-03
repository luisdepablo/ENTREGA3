%%%% SINTONIZACIÃ“N DE parametrosETROS Y OBTENCIÃ“N DE GRÃ?FICAS%%%%%
clear;
clc;
clf;
warning('off');

s=tf('s');

%parametrosetros del motor
K=2652.28;
p=64.986;
G=K/(s*(s+p));



%% Curvas de Mp
%zeta de 0.2 a 1
v=0.02;

Mp_des=1.115;
beta2=0.5;
betas=0:0.1:40;
zetas=0.2:0.2:1;
L=length(betas);

figure(1);
Legend=cell(length(zetas),1);
figure(1);
hold on;
n=0;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+(s/tau_d)+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_parametros(x,t,v);
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
betas=0:0.1:40;
zetas=1.2:0.2:2;
L=length(betas);
Legend=cell(length(zetas),1);
figure(2);
hold on;
n=0;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+(s/tau_d)+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_parametros(x,t,v);
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
betas=0:0.1:36;
zetas=0.2:0.2:2;
L=length(betas);

beta_xls=[];
n=0;

hold on;
for zeta=zetas
    tmp=[];
    for beta=betas
        [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
        H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^2*(s+p)+K*Kp*tau_d*(s^2+(s/tau_d)+1/(tau_d*tau_i)));
        [x,t]=step(H);
        [Mp,ts,tp,tr]=get_parametros(x,t,v);
        tmp=[tmp Mp];
        
    end
    
    n=n+1;
   [beta_inf,beta_sup]=get_betas(tmp,betas,Mp_des);
   beta_xls(n,1)=zeta;
   beta_xls(n,2)=beta_inf;
   beta_xls(n,3)=beta_sup;
   
   
    
end
filename='betas.csv';

writematrix(beta_xls,filename);

%% Calculo de ts
%para zetas entre 0.2 y 1

v=0.02;
load betas.csv
beta_xls=betas;

Legend=cell(6,1);
betas2=1:0.2:30;
zetas=0.2:0.2:1;
L=length(betas2);
ts_max=0.5;

n=0;
figure(3)

index=0;
index2=0;
hold on;


a=0;

for i=1:(length(zetas))
    
    n=n+1;
   
        tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,2);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;
            for beta2=betas2
                a=a+1;
%                 if tau_d==0
%                     pendiente=(tmp(end)-tmp(1))/(beta2-betas2(1));
%                     aux=pendiente*a+0.1;
%                     tmp=[tmp aux]
%                     
%                     continue;
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
                    
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_parametros(x,t,v);
                    tmp=[tmp ts];
                
              
            end
            
            
           
            
            
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;

           
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
        
            tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,3);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;         
            for beta2=betas2
                 [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
              
                H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
                [x,t]=step(H);
                [Mp,ts,tp,tr]=get_parametros(x,t,v);
                tmp=[tmp ts];
             end
            


%             beta2=betas2(find(ts==ts_max));
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;
            
            
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
    
end 
    
    
     

 
title('Curvas de t{s}')
plot(betas2,ones(L,1)*ts_max,'k--','LineWidth',1.2);
Legend{index+1}=strcat('t{s}');
legend(Legend,'FontSize',10);
xlabel('\beta_{2}');
ylabel('t{s}');
hold off;
saveas(gcf,'img/curvasTs1.fig')
saveas(gcf,'img/curvasTs1.png')

%% Calculo de ts
%para zetas entre 1.2 y 2
%debe ejecutarse justo despues del anterior para que se guarden las tuplas
v=0.02;



Legend=cell(11,1);
betas2=0.1:0.01:40;
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
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_parametros(x,t,v);
                    tmp=[tmp ts];
                    
              
            end
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;

           
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
            tmp=[];
        zeta=beta_xls(n,1);
        beta=beta_xls(n,3);
        
        if beta>0&&zeta>0
            index=index+1;
            index2=index2+1;
          
            for beta2=betas2
                
                    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
                    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
                    [x,t]=step(H);
                    [Mp,ts,tp,tr]=get_parametros(x,t,v);
                    tmp=[tmp ts];
             

            end
            
            [tupla]=get_tuplas(tmp,ts_max,betas2);
            tupla_xls(index2,1)=zeta;
            tupla_xls(index2,2)=beta;
            tupla_xls(index2,3)=tupla;
            
            Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(tupla));
            plot(betas2,tmp,'LineWidth',1.5);
        end
        
        
     
end
        filename2='tuplas.csv';         %zetas/betas/betas2
        writematrix(tupla_xls,filename2);  %0.2  /8.3  /18.6
                                        
 
title('Curvas de t{s}')
plot(betas2,ones(L,1)*ts_max,'k--','LineWidth',1.2);
Legend{index+1}=strcat('t_{s}max');
legend(Legend,'FontSize',10);
xlabel('\beta_{2}');
ylabel('t{s}');
hold off;
saveas(gcf,'img/curvasTs2.fig')
saveas(gcf,'img/curvasTs2.png')


%% respuesta al escalon
% 5 primeras
figure(5)
Legend=cell(5,1);
load tuplas.csv
load betas.csv
tupla_xls=tuplas;
beta_xls=betas;
t=0:0.001:4.0;
v=0.02;
hold on;

for i=1:5
    
zeta=tupla_xls(i,1);
beta=tupla_xls(i,2);
beta2=tupla_xls(i,3);

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    x=step(H,t);
    Legend{i}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1);
    
end

L=length(t);

title('Respuesta al escalón')
plot(t,ones(L,1)*(1+v),'k--','LineWidth',1.2);
plot(t,ones(L,1)*(1-v),'k--','LineWidth',1.2);
xl_ts=xline(0.5,'-','t_{s}Max');
xl_tr=xline(0.3,'-','t_{r}Max');
xl_ts.LabelOrientation='horizontal';
xl_tr.LabelOrientation='horizontal';
%plot(t,ones(L,1)*ts_max,'k--','LineWidth',1.2);

legend(Legend,'FontSize',10);
xlabel('t(s)');
ylabel('y(t)');
hold off;
saveas(gcf,'img/respuestas1.fig')
saveas(gcf,'img/respuestas1.png')



%% respuesta al escalon
%6 a 10 
figure(6)
Legend=cell(5,1);
hold on;
index=0;
t=0:0.001:1.5;
for i=6:10
    index=index+1;
    
zeta=tupla_xls(i,1);
beta=tupla_xls(i,2);
beta2=tupla_xls(i,3);

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    [x]=step(H,t);
    Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1);    
end

L=length(t);

title('Respuesta al escalón')
plot(t,ones(L,1)*(1+v),'k--','LineWidth',1);
plot(t,ones(L,1)*(1-v),'k--','LineWidth',1);
xl_ts=xline(0.5,'-','t_{s}Max');
xl_tr=xline(0.3,'-','t_{r}Max');
xl_ts.LabelOrientation='horizontal';
xl_tr.LabelOrientation='horizontal';
%plot(t,ones(L,1)*ts_max,'k--','LineWidth',1.2);

legend(Legend,'FontSize',10);
xlabel('t(s)');
ylabel('y(t)');

hold off;
saveas(gcf,'img/respuestas2.fig')
saveas(gcf,'img/respuestas2.png')




%% respuesta al escalon
%11 a 15 
figure(7)
Legend=cell(5,1);
t=0:0.001:1;
hold on;
index=0;
for i=11:15
    index=index+1;
    
zeta=tupla_xls(i,1);
beta=tupla_xls(i,2);
beta2=tupla_xls(i,3);

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    [x]=step(H,t);
    Legend{index}=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta_{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1);
    
end

L=length(t);

title('Respuesta al escalón')
plot(t,ones(L,1)*(1+v),'k--','LineWidth',1);
plot(t,ones(L,1)*(1-v),'k--','LineWidth',1);
xl_ts=xline(0.5,'-','t_{s}Max');
xl_tr=xline(0.3,'-','t_{r}Max');
xl_ts.LabelOrientation='horizontal';
xl_tr.LabelOrientation='horizontal';
%plot(t,ones(L,1)*ts_max,'k--','LineWidth',1.2);

legend(Legend,'FontSize',10);
xlabel('t(s)');
ylabel('y(t)');
hold off;
saveas(gcf,'img/respuestas2.fig')
saveas(gcf,'img/respuestas2.png')

%% respuesta al escalon
% cada grafica por separado

load tuplas.csv
tupla_xls=tuplas;

index=0;
for i=1:15
    index=index+1;
    figure(i);
zeta=tupla_xls(i,1);
beta=tupla_xls(i,2);
beta2=tupla_xls(i,3);

hold on;

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    [x,t]=step(H);
    Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
    plot(t,ones(length(t),1)*(1+v),'k--','LineWidth',1.2);
    plot(t,ones(length(t),1)*(1-v),'k--','LineWidth',1.2);
    xl_ts=xline(0.5,'-','t_{s}Max');
    xl_tr=xline(0.3,'-','t_{r}Max');
    yline(1.115,'-','M_{p}Max');
    xl_ts.LabelOrientation='horizontal';
    xl_tr.LabelOrientation='horizontal';
    
    legend(Legend);
hold off;
    
end

L=length(t);

%plot(t,ones(L,1)*ts_max,'k--','LineWidth',1.2);


title('Respuesta al escalón')
xlabel('t(s)');
ylabel('y(t)');


%% respuesta al escalon
% solo para tupla selecciona para punto 3
v=0.02;
load tuplas.csv
tupla_xls=tuplas;
t=0:0.001:1.5;


%tupla(n,i) n=2-> zeta=0.4 n=3-> zeta=0.6
   
    figure(1);
zeta=tupla_xls(2,1);
beta=tupla_xls(2,2);
beta2=tupla_xls(2,3)/4;
hold on;

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    [x]=step(H,t);
    Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
        plot(t,ones(length(t),1)*(1+v),'k--','LineWidth',1.2);
    plot(t,ones(length(t),1)*(1-v),'k--','LineWidth',1.2);
    xl_ts=xline(0.5,'-','t_{s}Max');
    xl_tr=xline(0.3,'-','t_{r}Max');
    yline(1.115,'-','M_{p}Max');
    xl_ts.LabelOrientation='horizontal';
    xl_tr.LabelOrientation='horizontal';
  

    
    legend(Legend);
    
    zeta=tupla_xls(2,1);
beta=tupla_xls(2,2);
beta2=tupla_xls(2,3);


    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+s/tau_d1+1/(tau_d1*tau_i)))/(s^3+K*Kp*tau_d*(s^2+s/tau_d+1/(tau_d*tau_i)));
    [x,t]=step(H);
    Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
        plot(t,ones(length(t),1)*(1+v),'k--','LineWidth',1.2);
    plot(t,ones(length(t),1)*(1-v),'k--','LineWidth',1.2);
    xl_ts=xline(0.5,'-','t_{s}Max');
    xl_tr=xline(0.3,'-','t_{r}Max');
    yline(1.115,'-','M_{p}Max');
    xl_ts.LabelOrientation='horizontal';
    xl_tr.LabelOrientation='horizontal';
  
hold off;
    
  

    





title('Respuesta al escalón')
xlabel('t(s)');
ylabel('y(t)');


