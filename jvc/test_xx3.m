clear all
i=0;
%thetab,phi0,psi0
thetab=pi/2;
phi0=pi/2;

%figure(2);clf;hold on;
figure(1);clf;hold on;
y=[];fy=[];x=[];fx=[];

%phi0=[pi/2,]
%for phi0=[pi/2,pi/2*0.9]
i=0;   
psi=-0.099; %0.8*(pi/2) / (pi/2);
Fy=0;
while psi<pi/2*0.9
   i=i+1;
   psi=psi+0.1;
   [Y(i),X(i),Fx(i),Fy(i),M(i)]=fxfym(thetab,phi0,psi);
   psi0(i)=psi;
end
figure(1);plot(X,Fy,'b',Y,Fy,'g',psi0/(pi/2),Fy,'r',psi0/(pi/2),3*psi0/(pi/2),'m');
%figure(1);plot(X,M,'b',Y,M,'g',psi0/(pi/2),M,'r');
xlabel('Y/L (green)    X/L (blue)    psi0/(pi/2) (red)    Linear (magenta)');ylabel('FyL^2/EI');title('nonlinear cantilever with both Fx and Fy');
%end
%figure(2);plot(X,Fx,'b',Y,Fx,'g');xlabel('Y/L (green)    X/L (blue)');ylabel('FxL^2/EI');title('nonlinear cantilever with both Fx and Fy');
%end
