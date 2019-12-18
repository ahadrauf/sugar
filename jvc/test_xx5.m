q_global=[0 0 0,0 1 0, 0 0 0, 0 1 0]';
Rp=[0 0 0]';
R=eye(3);
L=100e-6;
W=20e-6;
H=20e-6;

displaybeamold(q_global, Rp, R, L, W, H);

axis equal;
rotate3d on;
figure(1);
axis vis3d;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
shading interp;

xlabel('X');ylabel('Y');zlabel('Z');