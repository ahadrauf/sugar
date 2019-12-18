uses mumps.net

param flexure_length = 100u
param flexure_width = 2u
param flexure_length_2 = 20u
param flexure_width_2 = 2u
param mass_length = 50u
param mass_width = 50u
param F=100u

subnet crab_leg_part [a d] [N=*]
[
angle = pi/2*(2*floor(N/2)+1)
angle2 = angle + (2*N)*pi/2
anchor p1 [a] [l=10u w=10u oz=N*pi]

beam3d p1 [a b] [l=flexure_length w=flexure_width oz=N*pi] %fixed fleuxre
beam3dlink p1 [b c] [l=flexure_length_2 w=flexure_width_2 oz=-angle L2=mass_width/2-flexure_width_2/2 oz2=angle2] %flexure to mass

beam3dlink p1 [c d] [l=mass_length/2 w=mass_width oz=-angle L2=mass_width/2 oz2=angle] %proof mass

if N==0 | N==1
     f3d * [b] [F=F oz=pi/2]
]
f1 crab_leg_part p1 [w center] [N=0]
f2 crab_leg_part p1 [x center] [N=1]
f3 crab_leg_part p1 [y center] [N=2]
f4 crab_leg_part p1 [z center] [N=3]
% rigidlinkbeam p1 [w e] [l=100u w=10n h=10n oz=pi L1=0 L2=0]