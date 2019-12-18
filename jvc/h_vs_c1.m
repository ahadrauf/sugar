x=0;y=0;psi=0;X=0;Y=0;PSI=0;i=0;
for p=0.001:0.01/10:0.01
    i=i+1; [X(i),Y(i),PSI(i)]=hjelm_vs_cho_Fy_2(p); P(i)=p; 
    [x(i),y(i),psi(i)]=hjelmstad_nonlin_for2mat9(p);
    i
end
%figure(1);clf;plot(Y,P,'b',X,P,'r',PSI,P,'g');grid on;
%figure(2);clf;plot(y,P,'b',x,P,'r',psi,P,'g');grid on;
%figure(3);clf;plot(y-Y,P,'b',x-X,P,'r',psi-PSI,P,'g');grid on;


%figure(3);clf;plot(y(j),P(j),'b',x(j),P(j),'r',psi(j),P(j),'r');grid on;
%figure(1);clf;plot(Y(j),P(j),'b',X(j),P(j),'r',PSI(j),P(j),'g');grid on;
j=1:i;
I=((2e-6)^4)/12;E=165e9;L=20e-6;  
P=P*L*L/E/I;
figure(2);clf;plot(Y(j),P(j),'ob',X(j),P(j),'or',PSI(j),P(j),'og');grid on;hold on;
figure(2);plot(y(j),P(j),'b',-x(j),P(j),'r',psi(j),P(j),'g');grid on;

