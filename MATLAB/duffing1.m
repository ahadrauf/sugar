



%backbone
clear all; clc; clf;
warning('off');
format short;
global gamma q beta w0;
gamma = 0.5;
qq=[0.001:0.005:0.2]; %Forcing function
beta = 0.05;
w0=1;
%%%%%%%% initial conditions %%%%%%%%
a0=[0.1,0.1];
for iq=1:1:length(qq)
	q=qq(iq);
    options = optimset('TolFun',1e-15);
    ac =(fsolve(@(ac) funcAC(ac), a0, options));
    a0=[real(ac(1)),real(ac(2))];
    OM(iq,1)=((1+3*gamma/(16*w0^2)*ac(1)^2+9*gamma/(16*w0^2)*ac(1)^2)+sqrt(((1+3*gamma/(16*w0^2)*ac(1)^2+9*gamma/(16*w0^2)*ac(1)^2))^2-((1+3*gamma/(8*w0^2)*ac(1)^2)*(1+9*gamma/(8*w0^2)*ac(1)^2)+beta^2)));
    OM(iq,2)=((1+3*gamma/(16*w0^2)*ac(2)^2+9*gamma/(16*w0^2)*ac(2)^2)-sqrt(((1+3*gamma/(16*w0^2)*ac(2)^2+9*gamma/(16*w0^2)*ac(2)^2))^2-((1+3*gamma/(8*w0^2)*ac(2)^2)*(1+9*gamma/(8*w0^2)*ac(2)^2)+beta^2)));
    plot(real(OM(iq,2)), real(ac(2)),'ko');
    hold on;
    plot(real(OM(iq,1)), real(ac(1)),'k*');
end



