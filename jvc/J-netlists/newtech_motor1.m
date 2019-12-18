%p.qax=q(1);p.qay=q(2);p.qaz=q(3);   p.qaox=q(4);p.qaoy=q(5);p.qaoz=q(6);   p.F=20e-6;p.Foz=pi/2;p.Foy=0;  net=cho_load('newtech_motor1.m',p);q=cho_dc(net);figure(1);cho_display(net,q);
 
uses mumps.net

param F=0, Foz=0, Foy=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param LS_qx1=0, 	LS_qy1=0, 		LS_qz1=0
param LS_qox1=0, 	LS_qoy1=0, 	    LS_qoz1=0
param LS_qx2=0, 	LS_qy2=0, 		LS_qz2=0
param LS_qox2=0, 	LS_qoy2=0, 	    LS_qoz2=0
LS_L=100u
LS_W=2u
LS_H=2u

param RS_qx1=0, 	RS_qy1=0, 		RS_qz1=0
param RS_qox1=0, 	RS_qoy1=0, 	    RS_qoz1=0
param RS_qx2=0, 	RS_qy2=0, 		RS_qz2=0
param RS_qox2=0, 	RS_qoy2=0, 	    RS_qoz2=0
RS_L=100u
RS_W=2u
RS_H=2u

param RFS_qx1=0, 	RFS_qy1=0, 		RFS_qz1=0
param RFS_qox1=0, 	RFS_qoy1=0,	    RFS_qoz1=0
param RFS_qx2=0, 	RFS_qy2=0, 		RFS_qz2=0
param RFS_qox2=0, 	RFS_qoy2=0,	    RFS_qoz2=0
RFS_L=80u
RFS_W=2u
RFS_H=2u

param LFS_qx1=0, 	LFS_qy1=0, 		LFS_qz1=0
param LFS_qox1=0, 	LFS_qoy1=0,	    LFS_qoz1=0
param LFS_qx2=0, 	LFS_qy2=0,		LFS_qz2=0
param LFS_qox2=0, 	LFS_qoy2=0,	    LFS_qoz2=0
LFS_L=120u
LFS_W=2u
LFS_H=2u

param TFS_qx1=0, 	TFS_qy1=0, 		TFS_qz1=0
param TFS_qox1=0, 	TFS_qoy1=0,	    TFS_qoz1=0
param TFS_qx2=0, 	TFS_qy2=0, 		TFS_qz2=0
param TFS_qox2=0, 	TFS_qoy2=0,	    TFS_qoz2=0
TFS_L=100u
TFS_W=2u
TFS_H=2u

anchor p1 [AL]  [l=10u w=10u h=10u oz=pi/2] 
anchor p1 [AR]  [l=10u w=10u h=10u oz=pi/2] 
beam3d p1 [AL AM] [l=100u w=0.1u h=0.1u] 
anchor p1 [AM]  [l=10u w=10u h=10u oz=pi/2] 

LS  predisplacedbeam4 p1 [AL al][l=LS_L w=LS_W h=LS_H    ox=0 oy=0 oz=-pi/2   
    qx1  = LS_qx1 	    qy1  = LS_qy1 		qz1  = LS_qz1    
    qox1 = LS_qox1 	    qoy1 = LS_qoy1 	    qoz1 = LS_qoz1 
    qx2  = LS_qx2 	    qy2  = LS_qy2 		qz2  = LS_qz2    
    qox2 = LS_qox2     qoy2 = LS_qoy2	    qoz2 = LS_qoz2 ]

RS  predisplacedbeam4 p1 [AR ar][l=RS_L w=RS_W h=RS_H    ox=0 oy=0 oz=-pi/2   
    qx1  = RS_qx1 	    qy1  = RS_qy1 		qz1  = RS_qz1    
    qox1 = RS_qox1 	    qoy1 = RS_qoy1 	    qoz1 = RS_qoz1 
    qx2  = RS_qx2 	    qy2  = RS_qy2 		qz2  = RS_qz2    
    qox2 = RS_qox2     qoy2 = RS_qoy2	    qoz2 = RS_qoz2 ]

LFS predisplacedbeam4 p1 [al c][l=LFS_L w=LFS_W h=LFS_H    ox=0 oy=0 oz=0
    qx1  = LFS_qx1 	    qy1  = LFS_qy1 		qz1  = LFS_qz1    
    qox1 = LFS_qox1 	qoy1 = LFS_qoy1 	qoz1 = LFS_qoz1 
    qx2  = LFS_qx2 	    qy2  = LFS_qy2 		qz2  = LFS_qz2    
    qox2 = LFS_qox2     qoy2 = LFS_qoy2	    qoz2 = LFS_qoz2 ]

RFS predisplacedbeam4 p1 [c ar][l=RFS_L w=RFS_W h=RFS_H    ox=0 oy=0 oz=0
    qx1  = RFS_qx1 	    qy1  = RFS_qy1 		qz1  = RFS_qz1    
    qox1 = RFS_qox1 	qoy1 = RFS_qoy1 	qoz1 = RFS_qoz1 
    qx2  = RFS_qx2 	    qy2  = RFS_qy2 		qz2  = RFS_qz2    
    qox2 = RFS_qox2     qoy2 = RFS_qoy2	    qoz2 = RFS_qoz2 ]

TFS predisplacedbeam4 p1 [AM cc][l=TFS_L w=TFS_W h=TFS_H    ox=0 oy=0 oz=-(pi/2-pi/8)
    qx1  = TFS_qx1 	    qy1  = TFS_qy1 		qz1  = TFS_qz1    
    qox1 = TFS_qox1 	qoy1 = TFS_qoy1 	qoz1 = TFS_qoz1 
    qx2  = TFS_qx2 	    qy2  = TFS_qy2 		qz2  = TFS_qz2    
    qox2 = TFS_qox2     qoy2 = TFS_qoy2	    qoz2 = TFS_qoz2 ]

f3d * [cc][F=F oz=Foz oy=Foy]

%p.TFS_qy2=-19e-6+0*-1.75575703478355e-005;p.TFS_qx2=0.2e-6+0*-7.27173173882498e-006;p.TFS_qoz2=-0.285057767184866;p.F=0*-13e-6;net=cho_load('newtech_motor1.m',p);q=cho_dc(net);figure(2);cho_display(net,q);
%p.TFS_qy2=-19e-6+0*-1.75575703478355e-005;p.TFS_qx2=0.2e-6+0*-7.27173173882498e-006;p.TFS_qoz2=-0.285057767184866;p.F=-0e-6+0*-13e-6;net=cho_load('newtech_motor1.m',p);q=cho_dc(net);figure(2);cho_display(net,q);