i=0;
%2 bot
x=[0,1,2,3,4,5,6];
i=i+1;Displacement_y(i)=0.373915e-6;i=i+1;Displacement_y(i)=0.392426e-6;i=i+1;Displacement_y(i)=0.394318e-6;i=i+1;Displacement_y(i)=0.391793e-6;i=i+1;Displacement_y(i)=0.39075e-6;i=i+1;Displacement_y(i)=0.386464e-6;i=i+1;Displacement_y(i)=0.384877e-6;

%1 bot
x(8:14)=[0,1,2,3,4,5,6];
i=i+1;Displacement_y(i)=0.376014e-6;i=i+1;Displacement_y(i)=0.39596e-6;i=i+1;Displacement_y(i)=0.403277e-6;i=i+1;Displacement_y(i)=0.406073e-6;i=i+1;Displacement_y(i)=0.406912e-6;i=i+1;Displacement_y(i)=0.406907e-6;i=i+1;Displacement_y(i)=0.406541e-6;

%0.5 top
x(15:21)=[6,5,4,3,2,1,0];
i=i+1;Displacement_y(i)=0.415052e-6;i=i+1;Displacement_y(i)=0.414517e-6;i=i+1;Displacement_y(i)=0.413468e-6;i=i+1;Displacement_y(i)=0.41164e-6;i=i+1;Displacement_y(i)=0.407282e-6;i=i+1;Displacement_y(i)=0.396775e-6;i=i+1;Displacement_y(i)=0.376981e-6;

figure(1);
plot(x(1:7),Displacement_y(1:7),'*-b',x(8:14),Displacement_y(8:14),'*-g',x(15:21),Displacement_y(15:21),'*-r');
xlabel('Corner length [m],  blue = 2um mesh,   green = 1um mesh,   red = 0.5um mesh');
ylabel('Max Displacement [m],   F=30uN');
title('Deflection as a function of corner length. IntelliSuite simulation');

grid on;

L=20e-6;w=2e-6;h=2e-6;E=160e9;I=w^3*h/12;Fy=30e-6;
L3EI3=L^3/E/I/3;L2EI2=L^2/E/I/2;
L=[20,20+1/2,20+2/2,20+3/2,20+4/2,20+5/2,20+6/2]*1e-6;
D_y=Fy*L.^3/E/I/3; 
hold on;
plot(x(1:7),D_y,'o-b');


D_y=Fy*L.^3/E/I/3; 
L=((Displacement_y(15:21)/Fy*E*I*3).^(1/3) - 20e-6)
figure(2);
plot(x(15:21),L,'*-b');
grid on;
ylabel('Additional length that must be added to L. [m]');
xlabel('Corner length [m]');
title('Equivalent Roark beam. Length extension wrt corner length');