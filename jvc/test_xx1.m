clear all
i=0;
%thetab,phi0,psi0
thetab=pi/2;
phi0=pi/2;
thetab=pi/3;

figure(2);clf;hold on;
figure(1);clf;hold on;
y=[];fy=[];x=[];fx=[];
for thetab=0.1:.05:90*pi/180
   i=0;
   for psi=0.1:.05:90*pi/180
      i=i+1;
      [Y(i),X(i),Fx(i),Fy(i),M(i)]=fxfym(thetab,phi0,psi);
   end
%   figure(1);
%   plot3(X,Fy,Fx,'b',Y,Fy,Fx,'g');
figure(1);plot(X,Fy,'b',Y,Fy,'g');xlabel('Y/L (green)    X/L (blue)');ylabel('FyL^2/EI');title('nonlinear cantilever with both Fx and Fy');
figure(2);plot(X,Fx,'b',Y,Fx,'g');xlabel('Y/L (green)    X/L (blue)');ylabel('FxL^2/EI');title('nonlinear cantilever with both Fx and Fy');

x=[x;X];
fy=[fy;Fy];
fx=[fx;Fx];

end

figure(3);surf(x);
figure(4);surf(x,fy);
rotate3d on;

