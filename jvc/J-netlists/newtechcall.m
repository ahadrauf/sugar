clear all
F(1)=2e-6;
Foz(1)=pi/2;
F(2)=-10e-6;
Foz(2)=0;
F(3)=-2e-6;
Foz(3)=pi/2;
p.a=0;
for i=1:3
   
   p.F=F(i);
   p.Foz=Foz(i);
   
   net=cho_load('newtech3.m',p);
   q=cho_dc(net);
   figure(i);
   cho_display(net,q);
  
p.qax=q(lookup_coord(net,'a','x'));
p.qay=q(lookup_coord(net,'a','y'));
p.qaz=q(lookup_coord(net,'a','z'));
p.qaox=q(lookup_coord(net,'a','rx'));
p.qaoy=q(lookup_coord(net,'a','ry'));
p.qaoz=q(lookup_coord(net,'a','rz'));
p.qbx=q(lookup_coord(net,'b','x'));
p.qby=q(lookup_coord(net,'b','y'));
p.qbz=q(lookup_coord(net,'b','z'));
p.qbox=q(lookup_coord(net,'b','rx'));
p.qboy=q(lookup_coord(net,'b','ry'));
p.qboz=q(lookup_coord(net,'b','rz'));
p.qcx=q(lookup_coord(net,'c','x'));
p.qcy=q(lookup_coord(net,'c','y'));
p.qcz=q(lookup_coord(net,'c','z'));
p.qcox=q(lookup_coord(net,'c','rx'));
p.qcoy=q(lookup_coord(net,'c','ry'));
p.qcoz=q(lookup_coord(net,'c','rz'));

end

