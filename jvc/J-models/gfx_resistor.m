function gfx_resistor(tx,ty,tz,R,R1,R2,L,L1,L2,linewidth,nodewidth)
%jvclark - summer 2003
if 0
theta = pi/4;
s = sin(theta);
c = cos(theta);
R = [c -s 0; s c 0; 0 0 1];
tx = 0;ty = 0;tz = 0;
L = 100;L1=10;L2=50;
R1=R;R2=R;
%R=eye(3);
linewidth = 5;
color = 'k';
nodewidth = linewidth*3;
end

color = 'k';

%left link
q11=[0;0;0]; 
q12=R1*[L1;0;0]; 
q1=R*q11;
q2=R*q12;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%left lead
q21=q12; 
q22=q21+[L/5;0;0]; 
q1=R*q21;
q2=R*q22;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%resistor (saw wave)
xL=L/20;yL=L/5;
q31 = q22;
q32 = q31 + [xL;yL;0];
q1 = R*q31;
q2 = R*q32;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q41 = q32;
q42 = q41 + [2*xL;-2*yL;0];
q1 = R*q41;
q2 = R*q42;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q51 = q42;
q52 = q51 + [2*xL;2*yL;0];
q1 = R*q51;
q2 = R*q52;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q61 = q52;
q62 = q61 + [2*xL;2*-yL;0];
q1 = R*q61;
q2 = R*q62;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q71 = q62;
q72 = q71 + [2*xL;2*yL;0];
q1 = R*q71;
q2 = R*q72;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q81 = q72;
q82 = q81 + [2*xL;-2*yL;0];
q1 = R*q81;
q2 = R*q82;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q91 = q82;
q92 = q91 + [xL;yL;0];
q1 = R*q91;
q2 = R*q92;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%left lead
q101=q92; 
q102=q101+[L/5;0;0]; 
q1=R*q101;
q2=R*q102;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%left link
q111=q102; 
q112=q111+R2*[L2;0;0]; 
q1=R*q111;
q2=R*q112;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q2(1),ty+q2(2),tz+q2(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

