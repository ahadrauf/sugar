%fixfixnonlin6.m 

%running commands
%net=cho_load('fixfixlin6.m');dq=cho_dc(net);figure(1);cho_display(net,dq);dy=dq(lookup_coord(net,'b','y'))
%tic;net=cho_load('fixfixlin6.m');load=toc,tic;dq=cho_dc(net);dc=toc,figure(1);tic;cho_display(net,dq);dsp=toc,dof=net.dof,dy=dq(lookup_coord(net,'b','y'))
%dy=dq(lookup_coord(net,'b','y'))

%process file
uses mattmirrorprocess.m %E=185GPa
r=0
anchor p1 [a][l=10u w=10u h=10u oz=180+r]   
beam3d p1 [a b][l=100u w=6u h=6u oz=0+r]
beam3d p1 [b c][l=100u w=6u h=6u oz=0+r]
anchor p1 [c][l=10u w=10u h=10u oz=0+r]   

f3d * [b][F=2*1.062770294762911e-003 oz=90+r]

%[f,y]=fixfixFyY(50e-4)
%f = 1.062770294762911e-003
%y = 3.547462259020234e-006
%m = 3.501590893889557e-008

