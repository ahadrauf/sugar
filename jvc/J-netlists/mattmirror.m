%mattmirror.m     mirror/scanner design

%things to vary
momentarm=25u %25u, moment in x-direction in the hinge assembly by the ring
wbond=25u %25u, width of the bond beam
lbond=400u %40u, length of the bond beam
ltorsion=125u %25u, length of the ring torsions
wtorsion=10u %10u, width of the ring torsions
l2tether=35u %35u, thin beam length of the tether

%process file
uses mattmirrorprocess.m %E=185GPa

%constants
nfingers=9+2 %number of fingers
ltorsion2=25u %25u, length of the pull torsions
wtorsion2=20u %10u, width of the pull torsions
cmblink=200u %length of the link between combs
wthick=10u %minimum width of the thick beam
hbond=2u %2-4um thickness of the bond beam
hthick=50u %30-50um thickness of the thick beam
hthin=5u %3-5uthickness of the thin beams
l1tether=400u %thick beam length of the tether
s=0.6u %overlap offset
lcap=350u %length of capacitor links
lcapf=1400u %length of capacitor finger holders
crabparallel=40u %lenght of the parallel crab attachment
crabnormal=40u %lenght of the normal crab attachment
lcrab=400u %length of long crab part
lanchor=60u %anchor length
wanchor=60u %anchor width
hanchor=60u %anchor thickness
%ring constants
R=500u %radius of the ring
tr1=7.5 %angle increment of segment
tr2=15 %angle increment of segment
pi=3.141592653589793 %pi
A=R*sin(pi/2/3) %value for segment
B=R-R*cos(pi/2/3) %value for segment
Lring=sqrt(A*A+B*B)/2+1u %length of each ring segment

subnet base2hinge [d0 a2][M=* L=* W=* H=* OZ=* H2=* W2=* s=s] %a1 back mid, a2 mid-torsion
[ 
   s=0.2u
   beam3d parent [a1 b1][l=M   w=W   h=H  oz=OZ] %side1 extension
   beam3d parent [b1 b2][l=L   w=W   h=H  oz=OZ] %side1
   beam3d parent [b2 d0][l=L   w=W   h=H  oz=OZ-90] %back
   beam3d parent [d0 b3][l=L   w=W   h=H  oz=OZ-90] %back
   beam3d parent [b4 b3][l=L   w=W   h=H  oz=OZ] %side2
   beam3d parent [b5 b4][l=M   w=W   h=H  oz=OZ] %side2 extension   
   beam3d parent [b1 c1][l=H/2 w=W-s h=H  oy=90] %side1 torsion link dn
   beam3d parent [b4 c2][l=H/2 w=W-s h=H  oy=90] %side2 torsion link dn   
   beam3d parent [c1 a2][l=L   w=25u h=H2 oz=OZ-90] %torsion 1
   beam3d parent [a2 c2][l=L   w=25u h=H2 oz=OZ-90] %torsion 2
] 
subnet base3hinge1[a1 a2 a3][M=* L=* W=* H=* OZ=* H2=* W2=* s=s] %a1 side1 mid, a2 mid-torsion, a3 side2 top
[ 
   s=0.2u
   beam3d parent [a1 b1][l=M   w=W   h=H  oz=OZ] %side1 extension
   beam3d parent [b1 b2][l=25u w=W   h=H  oz=OZ] %side1
   beam3d parent [b2 b3][l=L*2 w=15u h=H  oz=OZ-90] %back
   beam3d parent [b4 b3][l=25u w=W   h=H  oz=OZ] %side2
   beam3d parent [b5 b4][l=M   w=W   h=H  oz=OZ] %side2 extension   
   beam3d parent [b1 c1][l=H/2 w=W-s h=H  oy=90] %side1 torsion link dn
   beam3d parent [b5 c2][l=H/2 w=W-s h=H  oy=90] %side2 torsion link dn   
   beam3d parent [b5 a3][l=H/2 w=W-s h=H  oy=-90] %side2 link up
   beam3d parent [c1 a2][l=L   w=W2  h=H2 oz=OZ-90] %torsion 1
   beam3d parent [a2 c2][l=L   w=W2  h=H2 oz=OZ-90] %torsion 2
] 
subnet tether [a1 a3][L1=* L2=* W=* H=* OZ=* H2=* s=s] %a1 side1 mid, a2 mid-torsion, a3 side2 top
[ 
   beam3d parent [b1 a1][l=H/2 w=W-s h=H/2 oy=90] %side1 torsion link dn
   beam3d parent [b1 b2][l=H/2 w=W-s h=H/2 oy=-90] %side2 torsion link up   
   beam3d parent [b1 b3][l=L1  w=W   h=H   oz=OZ] %side1 extension
   beam3d parent [b3 a2][l=H/2 w=W-s h=H/2 oy=90] %side1 torsion link dn
   beam3d parent [b3 b4][l=H/2 w=W-s h=H/2 oy=-90] %side2 torsion link up      
   beam3d parent [a2 a3][l=L2  w=W   h=H2  oz=OZ] %side1
] 
subnet comb [attach][platel=* platew=* fingerl=* fingerw=* nfinger=*]
[
   nfinger_half = (nfinger-1)/2
   platel_segment = platel/(nfinger-1)
   % Center finger
   beam3d parent [attach leg1] [l=fingerl w=fingerw]
   % First lower finger
   beam3d parent [attach lattach(1)] [l=platel_segment w=platew oz=-90]
   beam3d parent [lattach(1) lend(1)] [l=fingerl w=fingerw]
   % First upper finger
   beam3d parent [attach uattach(1)] [l=platel_segment w=platew oz=90]
   beam3d parent [uattach(1) uend(1)] [l=fingerl w=fingerw]
   for k = 2:nfinger_half-1 [
      % kth lower finger
      beam3d parent [lattach(k-1) lattach(k)] [l=platel_segment w=platew oz=-90]
      beam3d parent [lattach(k) lend(k)] [l=fingerl w=fingerw]
      % kth upper finge
      beam3d parent [uattach(k-1) uattach(k)] [l=platel_segment w=platew oz=90]
      beam3d parent [uattach(k) uend(k)] [l=fingerl w=fingerw]
   ]
]
subnet ring [r1 r2][L=* theta=* W=* H=* tr1=*] %a1 side1 mid, a2 mid-torsion, a3 side2 top
[ 
   beam3d parent [r1  r11][l=L w=W h=H oz=tr1+theta+tr2*0]
   beam3d parent [r11 r12][l=L w=W h=H oz=tr1+theta+tr2*1]
   beam3d parent [r12 r13][l=L w=W h=H oz=tr1+theta+tr2*2]
   beam3d parent [r13 r14][l=L w=W h=H oz=tr1+theta+tr2*3]
   beam3d parent [r14 r15][l=L w=W h=H oz=tr1+theta+tr2*4]
   beam3d parent [r15  r2][l=L w=W h=H oz=tr1+theta+tr2*5]
]
subnet mainring [r1 r3][]
[
   %ring
   beam3d p1 [r0 r1][l=R w=wthick h=hthick oz=-90]
   beam3d p1 [r0 r2][l=R w=wthick h=hthick oz=0]
   beam3d p1 [r0 r3][l=R w=wthick h=hthick oz=90]
   beam3d p1 [r0 r4][l=R w=wthick h=hthick oz=180]
   ring p1 [r1 r2][L=Lring W=wthick H=hthick tr1=tr1 theta=0]
   ring p1 [r2 r3][L=Lring W=wthick H=hthick tr1=tr1 theta=90]
   ring p1 [r3 r4][L=Lring W=wthick H=hthick tr1=tr1 theta=180]
   ring p1 [r4 r1][L=Lring W=wthick H=hthick tr1=tr1 theta=270]
]
subnet ringattachment [r1 cmb6][]
[
   %ring attachment
   base3hinge1 p1 [r1 t1 t2][L=ltorsion W=wthick H=hthick OZ=0 H2=hthin W2=wtorsion M=momentarm]
   tether p1 [t3 t1][L1=l1tether L2=l2tether W=wthick H=hthick OZ=0 H2=hthin]
   beam3d p1 [t2 A2][l=lbond w=wbond h=hbond oz=-90] %bond beam
   anchor p1 [A2][l=lanchor w=wanchor h=hanchor oy=90]
   base2hinge p1 [b1 t3][L=ltorsion2 W=wthick H=hthick OZ=180 H2=hthin W2=wtorsion2 M=momentarm]
   beam3d p1 [cmb1 b1][l=lcap w=wthick*2 h=hthick oz=0]
   %first crab
   beam3d p1 [cmb1 crab1][l=crabparallel w=wthick*2 h=hthick oz=-90]
   beam3d p1 [crab1 crab2][l=crabnormal w=wthick h=hthick oz=0]
   beam3d p1 [crab2 crab3][l=lcrab w=wthick h=hthick oz=-90]
   anchor p1 [crab3][l=lanchor w=wanchor h=hanchor]
   beam3d p1 [cmb1 crab4][l=crabparallel w=wthick*2 h=hthick oz=90]
   beam3d p1 [crab4 crab5][l=crabnormal w=wthick h=hthick oz=0]
   beam3d p1 [crab5 crab6][l=lcrab w=wthick h=hthick oz=90]
   anchor p1 [crab6][l=lanchor w=wanchor h=hanchor]
   %comb 
   comb p1   [cmb1][platew=wthick*2 platel=lcapf fingerl=100u fingerw=wthick nfinger=nfingers oz=180]
   beam3d p1 [cmb1 cmb2][l=cmblink w=wthick h=hthick oz=180]
   comb p1   [cmb2][platew=wthick*2 platel=lcapf fingerl=100u fingerw=wthick nfinger=nfingers oz=180]
   beam3d p1 [cmb2 cmb3][l=cmblink w=wthick h=hthick oz=180]
   comb p1   [cmb3][platew=wthick*2 platel=lcapf fingerl=100u fingerw=wthick nfinger=nfingers oz=180]
   beam3d p1 [cmb3 cmb4][l=cmblink w=wthick h=hthick oz=180]
   comb p1   [cmb4][platew=wthick*2 platel=lcapf fingerl=100u fingerw=wthick nfinger=nfingers oz=180]
   beam3d p1 [cmb4 cmb5][l=cmblink w=wthick h=hthick oz=180]
   comb p1   [cmb5][platew=wthick*2 platel=lcapf fingerl=100u fingerw=wthick nfinger=nfingers oz=180]
   beam3d p1 [cmb5 cmb6][l=cmblink w=wthick h=hthick oz=180]
   %end right crab
   beam3d p1 [cmb6 crab12][l=crabparallel w=wthick h=hthick oz=-90]
   beam3d p1 [crab12 crab22][l=crabnormal w=wthick h=hthick oz=180]
   beam3d p1 [crab22 crab32][l=lcrab w=wthick h=hthick oz=-90]
   anchor p1 [crab32][l=lanchor w=wanchor h=hanchor oz=180]
   beam3d p1 [cmb6 crab42][l=crabparallel w=wthick h=hthick oz=90]
   beam3d p1 [crab42 crab52][l=crabnormal w=wthick h=hthick oz=180]
   beam3d p1 [crab52 crab62][l=lcrab w=wthick h=hthick oz=90]
   anchor p1 [crab62][l=lanchor w=wanchor h=hanchor oz=180]   
]

mainring p1 [R L][] 
ringattachment p1 [R cmbR][] 
ringattachment p1 [L cmbL][oz=180] 

f3d * [cmbR][F=-100000u]
%f3d * [cmbL][F= 1000u]

