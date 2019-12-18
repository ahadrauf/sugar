% net=cho_load('newtech4a.m',p);q=cho_dc(net);figure(1);cho_display(net);
%p.F=40e-6;p.Foz=0;p.Foy=0;  p.qax=0;p.qay=100e-6;p.qaz=0;  p.qaox=0;p.qaoy=0;p.qaoz=0;  net=cho_load('newtech4a.m',p);q=cho_dc(net);figure(1);cho_display(net,q);
%p.F=-00e-6;p.Foz=0;p.Foy=0;  p.qax=0;p.qay=10e-6;p.qaz=0;  p.qbox=0;p.qboy=0;p.qboz=pi*1.5;  p.qaox=0;p.qaoy=0;p.qaoz=pi/4;  net=cho_load('newtech4a.m',p);q=cho_dc(net);figure(1);cho_display(net,q);
%p.F=40e-6;p.Foz=0;p.Foy=0;  net=cho_load('newtech4a.m',p);q=cho_dc(net);figure(1);cho_display(net,q);

%p.qax=q(1);p.qay=q(2);p.qaz=q(3);   p.qaox=q(4);p.qaoy=q(5);p.qaoz=q(6);   p.F=20e-6;p.Foz=pi/2;p.Foy=0;  net=cho_load('newtech4a.m',p);q=cho_dc(net);figure(1);cho_display(net,q);
 
uses mumps.net

param qax=0,qay=0,qaz=0,qaox=0,qaoy=0,qaoz=0
param qbx=0,qby=0,qbz=0,qbox=0,qboy=0,qboz=0
param qcx=0,qcy=0,qcz=0,qcox=0,qcoy=0,qcoz=0
param F=0, Foz=0, Foy=0

L=100u
W=2u
H=2u

anchor            p1 [A]  [l=10u w=10u h=10u oz=pi] 
ox=0
oy=0
oz=0

predisplacedbeam4 p1 [A a][l=L w=W h=H    ox=ox oy=oy oz=oz   
qox1=0    	qoy1=0    	qoz1=0      
qx2=qax 		qy2=qay 		qz2=qaz    
qox2=qaox 	qoy2=qaoy 	qoz2=qaoz ] 

%predisplacedbeam4 p1 [a b][l=L w=W h=H    ox=ox oy=oy oz=oz   qox1=qaox qoy1=qaoy qoz1=qaoz   qx2=qbx qy2=qby qz2=qbz    qox2=qbox qoy2=qboy qoz2=qboz ] 

%predisplacedbeam4 p1 [b c][l=L w=W h=H qox1=0.01*(qbox) qoy1=0.01*(qboy) qoz1=0.01*(qboz)   qx2=qcx qy2=qcy qz2=qcz qox2=qcox qoy2=qcoy qoz2=qcoz ] 

f3d * [a][F=F oz=Foz oy=Foy]

