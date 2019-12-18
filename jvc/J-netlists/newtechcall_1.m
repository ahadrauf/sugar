clear all
  F(1)=5e-6;
Foz(1)=pi/2;
  M(1)=0;

  F(2)=0;
Foz(2)=0;
  M(2)=-14e-10 ;


F(3)=1e-6;
Foz(3)=0;
  M(3)=0;

  F(4)=0;
Foz(4)=0;
  M(4)=0;

for i=1:4
   
   p.F=F(1);
   p.Foz=Foz(i);
   p.M=M(i);
   
   
   net=cho_load('newtech1.m',p);
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

