clear all
K=[  7, -3,   0;
     -3,  5,  -2;
     0,  -2,   2];
M=[ 3,0,0;
    0,2,0;
    0,0,1];

[e,ww]=eig(K,M);
D=2*0.02*sqrt(ww)*M;

F10=0.6;
%4.2336e-3*600;
F20=-4.42926e-3*600*0;
F30=3.3105e-3*600*0;
F=[F10;F20;F30];

I=eye(3,3);
O=zeros(3,3);
A=[-M\D , -M\K; I, O];
B=[M\F;0;0;0];
C=[1,0,0,0,0,0];
D=[0;0;0];


j=sqrt(-1);
i=0;
for w=0:0.01:6
    i=i+1;
    s=j*w;
    H=C*inv(s*[I,O;O,I]-A)*B;
    y0=abs(H)*F10;
    Y0(i)=y0(1);
    theta(i)=real(H)\imag(H);
    W(i)=w;
end
figure(1);
hold on;
plot(theta,Y0);
grid on;
