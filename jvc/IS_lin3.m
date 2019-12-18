L=100e-6;w=2e-6;h=2e-6;E=160e9;I=w^3*h/12;
L3EI3=L^3/E/I/3;L2EI2=L^2/E/I/2;

j=0;i=0;i=i+1;

j=j+1;T(i,j)=03.64;   meshing(i,j)=100e-6;  Fy(i,j)=10e-6;   Displacement_yt(i,j)=12.1205e-6;   Displacement_yb(i,j)=12.1203e-6;   Displacement_xt(i,j)=-0.235603e-6;   Displacement_xb(i,j)=4.94275e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
%j=j+1;T(i,j)=03.64;   meshing(i,j)=100e-6;  Fy(i,j)=10e-6;   Displacement_yt(i,j)=12.1205e-6;   Displacement_yb(i,j)=12.1203e-6;   Displacement_xt(i,j)=-0.235603e-6;   Displacement_xb(i,j)=4.94275e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));

j=j+1;T(i,j)=03.76;   meshing(i,j)=100e-6;  Fy(i,j)=20e-6;   Displacement_yt(i,j)=24.2411e-6;   Displacement_yb(i,j)=24.2406e-6;   Displacement_xt(i,j)=-0.471205e-6;   Displacement_xb(i,j)=9.88549e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=03.67;   meshing(i,j)=100e-6;  Fy(i,j)=30e-6;   Displacement_yt(i,j)=36.3616e-6;   Displacement_yb(i,j)=36.3609e-6;   Displacement_xt(i,j)=-0.706808e-6;   Displacement_xb(i,j)=14.8282e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.01;   meshing(i,j)=100e-6;  Fy(i,j)=40e-6;   Displacement_yt(i,j)=48.4821e-6;   Displacement_yb(i,j)=48.4812e-6;   Displacement_xt(i,j)=-0.942411e-6;   Displacement_xb(i,j)=19.771e-6;    D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.04;   meshing(i,j)=100e-6;  Fy(i,j)=50e-6;   Displacement_yt(i,j)=60.6026e-6;   Displacement_yb(i,j)=60.6015e-6;   Displacement_xt(i,j)=-1.17801e-6;    Displacement_xb(i,j)=24.7137e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=0;i=i+1;
j=j+1;T(i,j)=03.82;   meshing(i,j)=50e-6;   Fy(i,j)=10e-6;   Displacement_yt(i,j)=22.185e-6;    Displacement_yb(i,j)=22.1847e-6;   Displacement_xt(i,j)=-0.353291e-6;   Displacement_xb(i,j)=7.41511e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.22;   meshing(i,j)=50e-6;   Fy(i,j)=20e-6;   Displacement_yt(i,j)=29.5762e-6;   Displacement_yb(i,j)=29.5759e-6;   Displacement_xt(i,j)=-0.471005e-6;   Displacement_xb(i,j)=9.88685e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
Y_relative_error(i,j-1)=Y_relative_error(i,j);
j=j+1;T(i,j)=03.96;   meshing(i,j)=50e-6;   Fy(i,j)=30e-6;   Displacement_yt(i,j)=44.3643e-6;   Displacement_yb(i,j)=44.3639e-6;   Displacement_xt(i,j)=-0.706508e-6;   Displacement_xb(i,j)=14.8303e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.04;   meshing(i,j)=50e-6;   Fy(i,j)=40e-6;   Displacement_yt(i,j)=59.1524e-6;   Displacement_yb(i,j)=59.1518e-6;   Displacement_xt(i,j)=-0.942011e-6;   Displacement_xb(i,j)=19.7737e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.21;   meshing(i,j)=50e-6;   Fy(i,j)=50e-6;   Displacement_yt(i,j)=73.9405e-6;   Displacement_yb(i,j)=73.9398e-6;   Displacement_xt(i,j)=-1.17751e-6;    Displacement_xb(i,j)=24.7171e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=0;i=i+1;
j=j+1;T(i,j)=04.29;   meshing(i,j)=25e-6;   Fy(i,j)=10e-6;   Displacement_yt(i,j)=15.5252e-6;    Displacement_yb(i,j)=15.5251e-6;   Displacement_xt(i,j)=-0.235585e-6;  Displacement_xb(i,j)=4.94558e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.16;   meshing(i,j)=25e-6;   Fy(i,j)=20e-6;   Displacement_yt(i,j)=31.0505e-6;    Displacement_yb(i,j)=31.0501e-6;   Displacement_xt(i,j)=-0.47117e-6;   Displacement_xb(i,j)=9.89116e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.10;   meshing(i,j)=25e-6;   Fy(i,j)=30e-6;   Displacement_yt(i,j)=46.5757e-6;    Displacement_yb(i,j)=46.5752e-6;   Displacement_xt(i,j)=-0.706755e-6;  Displacement_xb(i,j)=14.8367e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.10;   meshing(i,j)=25e-6;   Fy(i,j)=40e-6;   Displacement_xt(i,j)=-0.94234e-6;   Displacement_xb(i,j)=19.7823e-6;   Displacement_yt(i,j)=62.101e-6;     Displacement_yb(i,j)=62.1003e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.01;   meshing(i,j)=25e-6;   Fy(i,j)=50e-6;   Displacement_xt(i,j)=-1.17793e-6;   Displacement_xb(i,j)=24.7279e-6;   Displacement_yt(i,j)=77.6262e-6;    Displacement_yb(i,j)=77.6254e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=0;i=i+1;
j=j+1;T(i,j)=04.51;   meshing(i,j)=12.5e-6;   Fy(i,j)=10e-6;   Displacement_xt(i,j)=-0.235741e-6;   Displacement_xb(i,j)=4.94921e-6;   Displacement_yt(i,j)=15.7212e-6;   Displacement_yb(i,j)=15.721e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.74;   meshing(i,j)=12.5e-6;   Fy(i,j)=20e-6;   Displacement_xt(i,j)=-0.471482e-6;   Displacement_xb(i,j)=9.89842e-6;   Displacement_yt(i,j)=31.4423e-6;   Displacement_yb(i,j)=31.442e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.50;   meshing(i,j)=12.5e-6;   Fy(i,j)=30e-6;   Displacement_xt(i,j)=-0.707224e-6;   Displacement_xb(i,j)=14.8476e-6;   Displacement_yt(i,j)=47.1635e-6;   Displacement_yb(i,j)=47.163e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.87;   meshing(i,j)=12.5e-6;   Fy(i,j)=40e-6;   Displacement_xt(i,j)=-0.942965e-6;   Displacement_xb(i,j)=19.7968e-6;   Displacement_yt(i,j)=62.8847e-6;   Displacement_yb(i,j)=62.884e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=04.55;   meshing(i,j)=12.5e-6;   Fy(i,j)=50e-6;   Displacement_xt(i,j)=-1.17871e-6;    Displacement_xb(i,j)=24.7461e-6;   Displacement_yt(i,j)=78.6058e-6;   Displacement_yb(i,j)=78.605e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=0;i=i+1;
j=j+1;T(i,j)=05.47;   meshing(i,j)=6.25e-6;   Fy(i,j)=10e-6;   Displacement_xt(i,j)=-0.235915e-6;   Displacement_xb(i,j)=4.95307e-6;   Displacement_yt(i,j)=15.7767e-6;   Displacement_yb(i,j)=15.7765e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=05.41;   meshing(i,j)=6.25e-6;   Fy(i,j)=20e-6;   Displacement_xt(i,j)=-0.471829e-6;   Displacement_xb(i,j)=9.90614e-6;   Displacement_yt(i,j)=31.5534e-6;   Displacement_yb(i,j)=31.553e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=05.46;   meshing(i,j)=6.25e-6;   Fy(i,j)=30e-6;   Displacement_xt(i,j)=-0.707744e-6;   Displacement_xb(i,j)=14.8592e-6;   Displacement_yt(i,j)=47.3301e-6;   Displacement_yb(i,j)=47.3295e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=05.49;   meshing(i,j)=6.25e-6;   Fy(i,j)=40e-6;   Displacement_xt(i,j)=-0.943659e-6;   Displacement_xb(i,j)=19.8123e-6;   Displacement_yt(i,j)=63.1068e-6;   Displacement_yb(i,j)=63.1061e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=05.66;   meshing(i,j)=6.25e-6;   Fy(i,j)=50e-6;   Displacement_xt(i,j)=-1.17957e-6;    Displacement_xb(i,j)=24.7653e-6;   Displacement_yt(i,j)=78.8834e-6;   Displacement_yb(i,j)=78.8826e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=0;i=i+1;
j=j+1;T(i,j)=09.01;   meshing(i,j)=3.125e-6;   Fy(i,j)=10e-6;   Displacement_xt(i,j)=-0.236397e-6;   Displacement_xb(i,j)=4.9633e-6;    Displacement_yt(i,j)=15.8302e-6;   Displacement_yb(i,j)=15.83e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=08.48;   meshing(i,j)=3.125e-6;   Fy(i,j)=20e-6;   Displacement_xt(i,j)=-0.472793e-6;   Displacement_xb(i,j)=9.9266e-6;    Displacement_yt(i,j)=31.6604e-6;   Displacement_yb(i,j)=31.66e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=08.20;   meshing(i,j)=3.125e-6;   Fy(i,j)=30e-6;   Displacement_xt(i,j)=-0.70919e-6;    Displacement_xb(i,j)=14.8899e-6;   Displacement_yt(i,j)=47.4906e-6;   Displacement_yb(i,j)=47.4901e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=08.51;   meshing(i,j)=3.125e-6;   Fy(i,j)=40e-6;   Displacement_xt(i,j)=-0.945587e-6;   Displacement_xb(i,j)=19.8532e-6;   Displacement_yt(i,j)=63.3208e-6;   Displacement_yb(i,j)=63.3201e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));
j=j+1;T(i,j)=08.58;   meshing(i,j)=3.125e-6;   Fy(i,j)=50e-6;   Displacement_xt(i,j)=-1.18198e-6;    Displacement_xb(i,j)==24.8165e-6;  Displacement_yt(i,j)=79.151e-6;    Displacement_yb(i,j)=79.1501e-6;   D_y(i,j)=Fy(i,j)*L3EI3; phi(i,j)=Fy(i,j)*L2EI2;   Y_relative_error(i,j) = ((D_y(i,j) - Displacement_yt(i,j))/D_y(i,j));


%meshing refinement vs relative error
figure(1);clf;
%plot3(L./meshing,Fy,abs(Y_relative_error));
surfl((L./meshing),Fy,(Y_relative_error));
xlabel('Meshing [# of discretizations]');
ylabel('Force [N]');
zlabel('Relative Error');
title('Small deflection analysis: Mesh vs Force vs RelativeError wrt Roark');
grid on;
rotate3d;

figure(2);
%plot3(L./meshing,Fy,abs(Y_relative_error));
surfl((L./meshing),Fy,T);
xlabel('Meshing [# of discretizations]');
ylabel('Force [N]');
zlabel('Time [s]');
title('Small deflection analysis: Mesh vs Force vs Time');
grid on;
rotate3d;

%meshing refinement vs relative error
figure(3);clf;
j=5;
plot(L./meshing(:,j),abs(Y_relative_error(:,j)));
xlabel('Meshing [# of discretizations]');
ylabel('Relative Error');
title('Small deflection analysis: Mesh vs RelativeError wrt Roark');
grid on;
rotate3d;
