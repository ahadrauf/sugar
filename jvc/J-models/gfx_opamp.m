function gfx_opamp(L,tx,ty,tz,R,linewidth)
theta = pi/4;
s = sin(theta);
c = cos(theta);
R = [c -s 0; s c 0; 0 0 1];
tx = 0;ty = 0;tz = 0;
L = 100;L1=10;L2=50;
R1=R;R2=R;
R=eye(3);
figure(1); clf; hold on;

linewidth = 5;
nodewidth = linewidth*3;
color = 'k';

%left leads
q11=[0;0;0]; 
q12=[L/4;0;0]; 
q1=R*q11;
q2=R*q12;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

q21=[0;-L/2;0]; 
q22=[L/4;-L/2;0]; 
q1=R*q21;
q2=R*q22;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%triangle
q31=[L/4;L/4;0]; 
q32=[L/4;-3*L/4;0]; 
q1=R*q31;
q2=R*q32;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

q41=[L/4;L/4;0]; 
q42=[3*L/4;-1*L/4;0]; 
q1=R*q41;
q2=R*q42;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

q41=[L/4;-3*L/4;0]; 
q42=[3*L/4;-1*L/4;0]; 
q1=R*q41;
q2=R*q42;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

q51=[3*L/4; -1*L/4;0]; 
q52=[L;     -1*L/4;0]; 
q1=R*q51;
q2=R*q52;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q2(1),ty+q2(2),tz+q2(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

