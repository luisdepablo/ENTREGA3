function[tau_d1,tau_d2,tau_d,tau_i,Kp]=trans_parametros(p,K,beta,beta2,zeta)
Kp=((p^2)*(2*beta+1/zeta^2))/((beta2^2)*K);
tau_d1=(beta2*(beta+2))/(p*(2*beta+1/zeta^2));
tau_d2=-p/(K*Kp);
tau_i=(beta2*(zeta^2)*(2*beta+1/zeta^2))/(beta*p);
tau_d=tau_d1+tau_d2;

end