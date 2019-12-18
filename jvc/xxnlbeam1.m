%You first choose a beme
%Y/L=1-sqrt(4*E*I/(P*L^2))*[E(k)-E(k,phi)]
%Where
%E(k)=complete elliptic integral of the second kind
%= integral(0 to pi/2)   sqrt(1-(k*sin(t))^2)*dt
%E(k,phi)=incomplete elliptic integral of the second kind
%= integral(0 to phi)   sqrt(1-(k*sin(t))^2)*dt
%k=sqrt((1+sin(theta))/2)
%P*L^2/(E*I)=(F(k)-F(k,phi))^2
%Phi=asin(1/(k*sqrt(2)));
%F(k)=complete elliptic integral of the first kind
%= integral(0 to pi/2)   1/sqrt(1-(k*sin(t))^2)*dt
%F(k,phi)=incomplete elliptic integral of the first kind
%= integral(0 to phi)   1/sqrt(1-(k*sin(t))^2)*dt
clear all;
resolution=80; %number of data steps
i=0; %data counter
theta1=0.001; %initial tip angle 
theta2=pi/2-0.01; %final tip angle
for theta=theta1:(theta2-theta1)/resolution:theta2
   i=i+1; %increment data counter
   k=sqrt((1+sin(theta))/2); %elliptic integral parameter
   phi=asin(1/(k*sqrt(2))); %elliptic integral limit
   integrand1st=strrep('1./sqrt(1-(k.*sin(t)).^2)','k',num2str(k)); %integrand of the 1st kind 
   integrand2nd=strrep('sqrt(1-(k.*sin(t)).^2)','k',num2str(k)); %integrand of the 2nd kind 
   %[Fc,Ec] = ellipke(k^2); %complete elliptic integrals of the first and second kinds
   %complete elliptic integral of the first kind 
   Fc=quad8(inline(integrand1st),0,pi/2); 
   %complete elliptic integral of the second kind 
   Ec=quad8(inline(integrand2nd),0,pi/2); 
   %incomplete elliptic integral of the first kind 
   Fi=quad8(inline(integrand1st),0,phi); 
   %incomplete elliptic integral of the second kind 
   Ei=quad8(inline(integrand2nd),0,phi); 
   %F*L^2/(E*I) is the nondimensional force, plotted on the y-axis 
   FLLEI(i)=(Fc-Fi)^2; 
   %Y/L, the nondimensional lateral displacement, plotted on the x-axis
   YL(i)=1-sqrt(4/FLLEI(i))*(Ec-Ei); 
   %X/L, the nondimensional beam shortening, plotted on the x-axis 
   XL(i)=1-sqrt(2*sin(theta)/FLLEI(i)); 
   %tipangle, the nondimensional tip angle, plotted on the x-axis
   tipangle(i)=theta/(pi/2); 
   %the linear elementary beam solution
   YLinearL(i)=3*YL(i); 
   %nonlinear stiffness
   K(i)=FLLEI(i)/YL(i);
   K2(i)=(Fc-Fi)^2/(1-2*(Ec-Ei)/(Fc-Fi));

end 
figure(1); plot(tipangle,FLLEI,'b',XL,FLLEI,'g',YL,FLLEI,'r',YL,YLinearL,'m',YL,K2,'k'); 
xlabel('Y/L (red),   X/L (green),   Theta/(pi/2) (blue),   YLinear/L (magenta)');
ylabel('force=Fy*L^2/(E*I),    (black)=stiffness vs Y/L');
title('Large deflection theory, L=length, E=Youngs, I=moment inertia about z-axis');


%K=L^3/(E*I)*(Fc-Fi)^2/(1-2*(Ec-Ei)/(Fc-Fi));
%Y=L*(1-2*(Ec-Ei)/(Fc-Fi));
%F=E*I/L^2*(Fc-Fi)^2;
