function gfx_vvt(tx,ty,tz,R,R1,R2,L,L1,L2)
RR = eye(3);

theta = pi/2;
s = sin(theta);
c = cos(theta);
R = [c -s 0; s c 0; 0 0 1];
tx = 0;ty = 0;tz = 0;
L = 100;L1=50;L2=50;
R1=R;R2=R';
%R=eye(3);
figure(1); clf; hold on;

linewidth = 5;
color = 'k';
nodewidth = linewidth*3;
os = [2.5*L1;0;0]; 

%left leads
q81=[-2.5*L1;0;0] + os; 
q82=[-1.8*L1;0;0] + os; 
q1=RR*q81;
q2=RR*q82;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

q91=[-2.5*L1;L;0] + os; 
q92=[-1.8*L1;L;0] + os; 
q1=RR*q91;
q2=RR*q92;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%left link
q11=[0;0;0]; 
q12=R1*[L1;0;0]; 
q1=RR*(R*q11 + os);
q2=RR*(R*q12 + os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%left lead
q21=q12; 
q22=q21+[L/4;0;0]; 
q1=RR*(R*q21 + os);
q2=RR*(R*q22 + os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%voltage (circle)
X = L/4; res=40;
q=[];
for x = -X:2*X/res:X;
   y = sqrt(-x^2 + X^2);
   qx = x + L/4 + q22(1);
   qy = y + q22(2);
   q = [q,RR*(R*[qx;qy;0]+os)];
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = -X:2*X/res:X;
   y = -sqrt(-x^2 + X^2);
   qx = x + L/4 + q22(1);
   qy = y + q22(2);
   q = [q,RR*(R*[qx;qy;0]+os)];
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 

%left lead
q31=q12+[3*L/4;0;0]; 
q32=q31+[L/4;0;0]; 
q1=RR*(R*q31+os);
q2=RR*(R*q32+os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%left link
q41=q32; 
q42=q41+R2*[L2;0;0]; 
q1=RR*(R*q41+os);
q2=RR*(R*q42+os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q2(1),ty+q2(2),tz+q2(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%minus
q51=q22 - [-2*L/16;-L/16;0]; 
q52=q51 + [0;-L/8;0]; ; 
q1=RR*(R*q51+os);
q2=RR*(R*q52+os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
%plus
q61=q31-[L/16;0;0]; 
q62=q61-[L/8;0;0]; 
q1=RR*(R*q61+os);
q2=RR*(R*q62+os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
q61=q31+[-2*L/16;-L/16;0]; 
q62=q61-[0;-L/8;0]; 
q1=RR*(R*q61+os);
q2=RR*(R*q62+os);
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

