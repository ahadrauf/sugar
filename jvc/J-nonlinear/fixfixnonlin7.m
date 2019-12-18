%fixfixnonlin5.m 

%running commands
%format compact; format long e;
%net=cho_load('fixfixnonlin7.m');dq=cho_dc(net);figure(1);cho_display(net,dq);X=dq(lookup_coord(net,'b','x')),Y=dq(lookup_coord(net,'b','y')),RZ=dq(lookup_coord(net,'b','rz')) 
%net=cho_load('fixfixnonlin5.m');dq=cho_dc(net);figure(1);cho_display(net,dq);dy=dq(lookup_coord(net,'b','y'))
%tic;net=cho_load('fixfixnonlin5.m');load=toc,tic;dq=cho_dc(net);dc=toc,figure(1);tic;cho_display(net,dq);dsp=toc,dof=net.dof,dy=dq(lookup_coord(net,'b','y'))
%dy=dq(lookup_coord(net,'b','y'))

%process file
uses mattmirrorprocess.m %E=185GPa
r=0
r2=0
anchor    p1 [aa]   [l=10u w=10u h=10u oz=180+r]   
beam3dnl4 p1 [aa a] [l=50u w=6u  h=6u  oz=r2+r]
beam3dnl4 p1 [a b]  [l=50u w=6u  h=6u  oz=0+r]
beam3dnl4 p1 [b c]  [l=50u w=6u  h=6u  oz=0+r]
beam3dnl4 p1 [c cc] [l=50u w=6u  h=6u  oz=-r2+r]
anchor    p1 [cc]   [l=10u w=10u h=10u oz=0+r]   

f3d * [b][F=2*1.062770294762911e-003 oy=-90]
f3d * [b][F=2*1.062770294762911e-003 oz=90+r]

%[f,y]=fixfixFyY(300e-4)
%f = 5.223913491604605e-003
%y = 8.790989065775853e-006

%[f,y]=fixfixFyY(100e-4)
%f = 1.806269814656536e-003
%y = 5.030236143960421e-006

%[f,y]=fixfixFyY(50e-4)
%f = 1.062770294762911e-003
%y = 3.547462259020234e-006
%m = 3.501590893889557e-008

