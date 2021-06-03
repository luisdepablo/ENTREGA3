clc;clear;clf
K=2652.28;
p=64.986;
beta=35.24;
beta2=19.14;
zeta=0.2;
T=5e-3;
v=0.02;
s=tf('s');

%%


%%
[tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta);

%% datos del telelabo
clf;
figure(1)


file="test-MOTOR3REF";
Tr=readtable(file);
ar=table2array(Tr);
Legend=cell(3,1);

hold on;
t=0:0.01:5.01;

L=length(t);

plot(ar(:,1),ar(:,2));
Legend{1}=strcat('Referencia');

file="test-MOTOR3POS";
Tr=readtable(file);
ar=table2array(Tr);


plot(ar(:,1),ar(:,2));
Legend{2}=strcat('Real');







% simulacion


    H=(K*Kp*tau_d1*(s^2+(s/tau_d1)+1/(tau_d1*tau_i)))/(s^3+(p+K*Kp*tau_d)*s^2+(K*Kp)*s+(K*Kp/tau_i));
    x=step((1/s)*H,t);
    %Legend=strcat('\zeta=',num2str(zeta),', \beta=',num2str(beta),' \beta{2}=',num2str(beta2));
    plot(t,x,'LineWidth',1.5);
    
   
    xlim([0,4]);
    Legend{3}=strcat('Simulado');
  
hold off;
    title("Respuesta a la rampa");
    xlabel("tiempo (s)")
    ylabel("posicion (rad)")
    legend(Legend);
    saveas(gcf,'respuestaMotorRealRampa.fig')
%% grafica del error

% figure(2);
% 
% file="test-MOTOR3ERR";
% 
% Tr=readtable(file);
% ar=table2array(Tr);
% Legend=cell(2,1);
% hold on;
% plot(ar(:,1),ar(:,2));
% Legend{1}=strcat('error');
% 
% file="POS-REF1";
% Tr=readtable(file);
% ar=table2array(Tr);
% 
% 
% plot(ar(:,1),ar(:,2));
% Legend{2}=strcat('Real');
% 
% hold off;
% xlabel("tiempo (s)")
% ylabel("posicion (rad)")
% xlim([0,5]);
% legend(Legend);

%% ERROR
figure(2)
file="test-MOTOR3ERR";
Tr=readtable(file);
ar=table2array(Tr);

plot(ar(:,1),ar(:,2));
xlim([0,4]);

Legend=strcat('Error');
xlabel("tiempo (s)")
ylabel("posicion (rad)")
title("Error de la rampa");
legend(Legend);

saveas(gcf,'errorRealRampa.fig')