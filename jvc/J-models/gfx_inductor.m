function gfx_inductor(tx,ty,tz,R,R1,R2,L,L1,L2,linewidth,nodewidth)
%jvclark - summer 2003
test = 0;
if test
    theta = pi/4;
    s = sin(theta);
    c = cos(theta);
    R = [c -s 0; s c 0; 0 0 1];
    tx = 0;ty = 0;tz = 0;
    L = 100;L1=10;L2=50;
    R1=R;R2=R;
    %R=eye(3);
    linewidth = 2;
    nodewidth = 4;
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
q22=q21+[2*L/16;0;0]; 
q1=R*q21;
q2=R*q22;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

X = L/4; res=200;

%left circle
q=[];
a=4; b=-L/5;
q22 = [3*L/16;0;0] + q12;
for x = -X:X/res:X;
   y = sqrt(-a*x^2 + X^2); qx = x + 2*L/32 + q22(1); qy = y + q22(2);
   if y == real(y)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = 0:X/res:X;
   y = -sqrt(-a*x^2 + X^2); qx = x + 2*L/32 + q22(1); qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 


%middle (circle)
q=[];
for x = -X:2*X/res:X;
   y = sqrt(-a*x^2 + X^2);
   qx = x + 7*L/32 + q22(1);
   qy = y + q22(2);
   if y == real(y)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = -X:X/res:0;
   y = -sqrt(-a*x^2 + X^2); qx = x + 7*L/32 + q22(1); qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = 0:X/res:X;
   y = -sqrt(-a*x^2 + X^2); qx = x + 7*L/32 + q22(1); qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 


%right circle
q=[];
for x = -X:2*X/res:X;
   y = sqrt(-a*x^2 + X^2);
   qx = x + 12*L/32 + q22(1);
   qy = y + q22(2);
   if y == real(y)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = 0:X/res:X;
   y = -sqrt(-a*x^2 + X^2);qx = x + 12*L/32 + q22(1);qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = -X:X/res:0;
   y = -sqrt(-a*x^2 + X^2);qx = x + 12*L/32 + q22(1);qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 

%far right circle
q=[];
for x = -X:2*X/res:X;
   y = sqrt(-a*x^2 + X^2);
   qx = x + 17*L/32 + q22(1);
   qy = y + q22(2);
   if y == real(y)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 
q=[];
for x = -X:2*X/res:0
   y = -sqrt(-a*x^2 + X^2);qx = x + 17*L/32 + q22(1);qy = y + q22(2);
   if (y == real(y))&(y>b)
       q = [q,R*[qx;qy;0]];
   end
end
plot3(tx+q(1,:),ty+q(2,:),tz+q(3,:),color,'LineWidth',linewidth); 

%left lead
q31=q12+[27*L/32;0;0]; 
q32=q12+[L;0;0]; 
q1=R*q31;
q2=R*q32;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%left link
q41=q32; 
q42=q41+R2*[L2;0;0]; 
q1=R*q41;
q2=R*q42;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q2(1),ty+q2(2),tz+q2(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

