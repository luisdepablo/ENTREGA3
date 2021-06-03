clc;clear;clf
K=2652.28/23;
p=64.986;
beta=17.22;
beta2=16.42;
zeta=0.4;
T=5e-3;
v=0.02;
s=tf('s');



%% datos del telelabo
clf;
figure(1)

file="POS-REF1";
Tr=readtable(file);
ar=table2array(Tr);
Legend=cell(2,1);

hold on;
t=0:0.01:5.01;

L=length(t);

plot(ar(:,1),ar(:,2));
Legend{1}=strcat('Real');

%%simulacion

    [tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);
    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
    x=step(H,t);
    %Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
    plot(t,ones(L,1)*(1+v),'k--','LineWidth',1.2);
    plot(t,ones(L,1)*(1-v),'k--','LineWidth',1.2);
    xl_ts=xline(0.5,'-','t_{s}Max');
    xl_tr=xline(0.3,'-','t_{r}Max');
    yline(1.115,'-','M_{p}Max');
    xl_ts.LabelOrientation='horizontal';
    xl_tr.LabelOrientation='horizontal';
    xl_tr.LabelVerticalAlignment='middle';
    xlim([0,5]);
    Legend{2}=strcat('Simulado');
  
hold off;

    xlabel("tiempo (s)")
    ylabel("posicion (rad)")
    title("Respuesta al escalón");
    legend(Legend);
saveas(gcf,'respuestaMotorRealEscalon.fig')
%% grafica del error

figure(2);

file="test-MOTOR3ERR";

Tr=readtable(file);
ar=table2array(Tr);
Legend=cell(2,1);
hold on;
plot(ar(:,1),ar(:,2));
Legend{1}=strcat('Error');

file="POS-REF1";
Tr=readtable(file);
ar=table2array(Tr);


plot(ar(:,1),ar(:,2));
Legend{2}=strcat('Real');

hold off;
xlabel("tiempo (s)")
ylabel("posicion (rad)")
xlim([0,5]);
title("Respuesta del error");
legend(Legend);
saveas(gcf,'respuestaErrorMotorRealEscalon.fig')


%% grafica de respuesta mejorada
figure(3)
file="test-MOTOR3POS";
Tr=readtable(file);
ar=table2array(Tr)/3.14;
Legend=cell(2,1);
beta2=4.115;

[tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);


hold on;
t=0:0.00001:1;

L=length(t);

plot(ar(:,1),ar(:,2),'LineWidth',1.5);
Legend{1}=strcat('Real');


%%simulacion
    
    beta2=beta2/4;
    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
    x=step(H,t);
    %Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
    plot(t,ones(L,1)*(1+v),'k--','LineWidth',1.2);
    plot(t,ones(L,1)*(1-v),'k--','LineWidth',1.2);
    xl_ts=xline(0.5,'-','t_{s}Max');
    xl_tr=xline(0.3,'-','t_{r}Max');
    yline(1.115,'-','M_{p}Max');
    xl_ts.LabelOrientation='horizontal';
    xl_tr.LabelOrientation='horizontal';
    xlim([0,0.6]);
    Legend{2}=strcat('Simulado');
  
hold off;

    xlabel("tiempo (s)")
    ylabel("posicion (rad)")
    title("Respuesta al escalón");
    legend(Legend);
    
    saveas(gcf,'respuestaMotorRealEscalonMejorado.fig')
