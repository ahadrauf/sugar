L=100e-6;w=2e-6;h=2e-6;E=160e9;I=w^3*h/12;
L3EI3=L^3/E/I/3;L2EI2=L^2/E/I/2;

i=0;
i=i+1;
Fy(i)=50e-6; meshing(i)=100e-6;  Displacement_y(i)=57.1007e-6; D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);
i=i+1;
Fy(i)=50e-6; meshing(i)=50e-6;   Displacement_y(i)=70.6209e-6; D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);
i=i+1;
Fy(i)=50e-6; meshing(i)=25e-6;   Displacement_y(i)=75.4486e-6; D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);
i=i+1;
Fy(i)=50e-6; meshing(i)=12.5e-6; Displacement_y(i)=77.0801e-6; D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);
i=i+1;
Fy(i)=50e-6; meshing(i)=6.25e-6; Displacement_y(i)=77.6823e-6; D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);
i=i+1;
Fy(i)=50e-6; meshing(i)=3.125e-6;Displacement_y(i)=77.9153e-6;D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_y(i))/D_y);

i=i+1;
meshing(i)=100e-6;  Fy(i)=50e-6;   Displacement_yt(i)=60.6026e-6;   D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);
i=i+1;
meshing(i)=50e-6;   Fy(i)=50e-6;   Displacement_yt(i)=73.9405e-6;   D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);
i=i+1;
meshing(i)=25e-6;   Fy(i)=50e-6;   Displacement_yt(i)=77.6262e-6;   D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);
i=i+1;
meshing(i)=12.5e-6; Fy(i)=50e-6;   Displacement_yt(i)=78.6058e-6;   D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);
i=i+1;
meshing(i)=6.25e-6; Fy(i)=50e-6;   Displacement_yt(i)=78.8834e-6;   D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);
i=i+1;
meshing(i)=3.125e-6;Fy(i)=50e-6;   Displacement_yt(i)=79.151e-6;    D_y=Fy(i)*L3EI3; Y_relative_error(i) = ((D_y - Displacement_yt(i))/D_y);


L=100e-6;

%meshing refinement vs relative error
figure(1);clf;
j=1:6;
plot((L./meshing(j)),(Y_relative_error(j)),'-ob');
xlabel('Meshing [# of discretizations],  blue=no anchor, red=anchor');
%ylabel('Force [N]');
ylabel('RelativeError = (Roark-IntelliSuite)/Roark');
title('Cantilever w&w/o anchor. Linear analysis. Mesh vs RelativeError');
grid on;
rotate3d;

hold on
%figure(1);clf;
j=7:12;
plot((L./meshing(j)),(Y_relative_error(j)),'o-r');
xlabel('Meshing [# of discretizations],  blue=no anchor, red=anchor');
%ylabel('Force [N]');
ylabel('RelativeError = (Roark-IntelliSuite)/Roark');
title('Cantilever w&w/o anchor. Linear analysis. Mesh vs RelativeError');
grid on;
rotate3d;


figure(2);clf;
j=2:6;
plot((L./meshing(j)),(Displacement_y(j)-Displacement_y(j-1)),'o-b');
xlabel('Meshing [# of discretizations],  blue=no anchor, red=anchor');
%ylabel('Force [N]');
ylabel('Convergence = y(i)-y(i-1)  [m]');
title('Cantilever w&w/o anchor. Linear analysis. Mesh vs y(i)-y(i-1)');
grid on;
rotate3d;

hold on
j=8:12;
plot((L./meshing(j)),(Displacement_yt(j)-Displacement_yt(j-1)),'o-r');
xlabel('Meshing [# of discretizations],  blue=no anchor, red=anchor');
%ylabel('Force [N]');
ylabel('Convergence = y(i)-y(i-1)  [m]');
title('Cantilever w&w/o anchor. Linear analysis. Mesh vs y(i)-y(i-1)');
grid on;
rotate3d;

