function gfx_impedance(L,tx,ty,tz,R,linewidth,R1,R2,L1,L2) 
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
color = 'k';
nodewidth = linewidth*3;

%left dot link 
q11=[0;0;0]; 
q12=R1*[L1;0;0]; 
q1=R*q11;
q2=R*q12;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%leads
q21=q12;
q22=q21+[L/8;0;0];
q1=R*q21;
q2=R*q22;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

q31=q12+[7*L/8;0;0];
q32=q12+[L;0;0];
q1=R*q31;
q2=R*q32;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);

%box
q51=q22 + [0;L/16;0];
q52=q22 + [0;-L/16;0];
q1=R*q51;
q2=R*q52;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);
q61=q22 + [3*L/4;L/16;0];
q62=q22 + [3*L/4;-L/16;0];
q1=R*q61;
q2=R*q62;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);
q71=q51;
q72=q61;
q1=R*q71;
q2=R*q72;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);
q71=q52;
q72=q62;
q1=R*q71;
q2=R*q72;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);

%right dot link
q41=q32;
q42=q41+R2*[L2;0;0];
q1=R*q41;
q2=R*q42;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);
plot3(tx+q2(1),ty+q2(2),tz+q2(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

