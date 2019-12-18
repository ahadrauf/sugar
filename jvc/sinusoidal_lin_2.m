clear all;
t=0;
E=0.02;
w1=18.85;
w2=42.39;
w3=58.718;
w12=w1^2;
w22=w2^2;
w32=w3^2;
F10=4.2336e-3*600;
F20=-4.42926e-3*600;
F30=3.3105e-3*600;
F1=F10*sin(w1*t);
F2=F20*sin(w1*t);
F3=F30*sin(w1*t);

M=eye(3,3);
C=2*E*[w1,0,0;0,w2,0;0,0,w3];
K=[w12,0,0;0,w22,0;0,0,w32];
F=[F10,0,0;0,F20,0;0,0,F30];
I=eye(3,3);
O=zeros(3,3);
A=[C,K;I,O];
B=[F,O;O,O];
C=[I,O;O,O];
w=18;
s=sqrt(-1)*w;
H=C*inv(s*[I,O;O,I]-A)*B;

j=sqrt(-1);
i=0;
for w=0:80
    i=i+1;
    s=j*w;
    H=[s^2+2*E*w1*s+w12; 
       s^2+2*E*w2*s+w22; 
       s^2+2*E*w3*s+w32];
    y0=abs(H);
    Y0(i)=y0(3);
    theta(i)=real(H)\imag(H);
    W(i)=w;
end
figure(1);
hold on;
plot(W,Y0);
grid on;
