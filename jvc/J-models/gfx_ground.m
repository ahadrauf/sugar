function gfx_ground(tx,ty,tz,R,R1,R2,L,L1,L2,linewidth,nodewidth)
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

%left dot link 
q11=[0;0;0]; 
q12=R1*[L1;0;0]; 
q1=R*q11;
q2=R*q12;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 
plot3(tx+q1(1),ty+q1(2),tz+q1(3),'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',nodewidth); 

%leads
q21=q12;
q22=q21+[2.5*L/6;0;0];
q1=R*q21;
q2=R*q22;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth); 

%gap
q1=q22+[0;L/4;0];
q2=q22+[0;-L/4;0];
q1=R*q1;
q2=R*q2;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);

q1=q22+[L/16;L/7;0];
q2=q22+[L/16;-L/7;0];
q1=R*q1;
q2=R*q2;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);

q1=q22+[2*L/16;L/16;0];
q2=q22+[2*L/16;-L/16;0];
q1=R*q1;
q2=R*q2;
plot3(tx+[q1(1),q2(1)],ty+[q1(2),q2(2)],tz+[q1(3),q2(3)],color,'LineWidth',linewidth);