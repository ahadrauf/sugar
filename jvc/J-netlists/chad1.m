uses mumps.net

a=15u
b=50u
l=100u
w=200u
h=1.25u
x=100u
y=20u

subnet corner [d g][]
[
%upper left

anchor p1 [a][ l=a w=b h=h oz=-90]
anchor p1 [f][ l=a w=b h=h oz=180]

%cantilever
beam3d p1 [a b][l=x w=y h=h oz=-90 ]
beam3d p1 [e f][l=x w=y h=h oz=180 ]

beam3d p1 [b c][l=l w=w h=h oz=-90 ]
beam3d p1 [c d][l=l w=w h=h oz=-90 ]
beam3d p1 [c e][l=l w=w h=h oz=180 ]
beam3d p1 [c g][l=l w=w h=h oz=0 ]
]

subnet corner2 [d g][]
[
%upper left

anchor p1 [a][ l=a w=b h=h oz=-90]
anchor p1 [f][ l=a w=b h=h oz=0]

%cantilever
beam3d p1 [a b][l=x w=y h=h oz=-90 ]
beam3d p1 [e f][l=x w=y h=h oz=180 ]

beam3d p1 [b c][l=l w=w h=h oz=-90 ]
beam3d p1 [c d][l=l w=w h=h oz=-90 ]
beam3d p1 [c e][l=l w=w h=h oz=180 ]
beam3d p1 [c g][l=l w=w h=h oz=0 ]
]


subnet sidebeam [a c][]
[
beam3d p1 [a b][l=150u w=w h=h oz=0 ]
beam3d p1 [b c][l=150u w=w h=h oz=0 ]
beam3d p1 [b d][l=100u w=300u h=h oz=90]
]

subnet sidebeam2 [a c][]
[
beam3d p1 [a b][l=150u w=w h=h oz=0 ]
beam3d p1 [b c][l=150u w=w h=h oz=0 ]
beam3d p1 [b d][l=100u w=300u h=h oz=-90]
]



corner p1 [a b][]
sidebeam p1 [b c][]
sidebeam p1 [c d][]
sidebeam2 p1 [a f][oz=-90]
sidebeam2 p1 [f g][oz=-90]
corner2 p1 [d e][oz=-90]
sidebeam p1 [e h][oz=-90]
sidebeam p1 [h i][oz=-90]
corner p1 [i j][oz=180]
sidebeam p1 [j k][oz=180]
sidebeam p1 [k l][oz=180]
corner2 p1 [l g][oz=90]

beam3d p1 [c k][ l=800u w=600u h=h oz=-90]

beam3d p1 [f h][ l=800u w=600u h=h oz=0]


q=0.5u

f3d * [a][F=q oy=90]
f3d * [b][F=q oy=90]
f3d * [c][F=q oy=90]
f3d * [d][F=q oy=90]
f3d * [f][F=q oy=90]
f3d * [g][F=q oy=90]
f3d * [h][F=q oy=90]
f3d * [i][F=q oy=90]
f3d * [j][F=q oy=90]
f3d * [k][F=q oy=90]
f3d * [l][F=q oy=90]
f3d * [e][F=q oy=90]
