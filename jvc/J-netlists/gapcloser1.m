%net=cho_load('gapcloser1.m');figure(1);cho_display(net);


uses mumps.net
rot=pi/8
dw=10u
subnet gapcloser [a(1,1)][l=* w=* oz=* h=* gap=*]
[  
   for j=1:3
   [
      for i=1+2:5+2
      [
         anchor p1 [a(i,j)][l=40u w=40u h=h*4 oz=pi]
         beam2d p1 [a(i,j) b(i,j)][l=l w=5u oz=(2-j+1)*rot h=h]
         gap2d p1 [b(i,j) c(i,j) d(i,j) e(i,j)][gap=gap l=150u w1=w+(i*dw) w2=20u oz=oz h=h oz=(2-j+1)*rot]
         beam2d p1 [a(i,j) a(i+1,j)][oz=-pi/2 l=100u+w+(i*dw) w=0u h=0u]      
         
%         anchor p1 [d(i,j)][l=20u w=40u h=h*4 oz=-pi/2]
%         anchor p1 [e(i,j)][l=20u w=40u h=h*4 oz=-pi/2]
         
         anchor p1 [f(i,j)][l=40u w=40u h=h*4 oz=pi]
         beam2d p1 [f(i,j) d(i,j)][l=l+((gap+(w+(i*dw))/2+20u/2)*sin((2-j+1)*rot))/cos((2-j+1)*rot) w=2u oz=(2-j+1)*rot h=h/2]
      ]
      beam2d p1 [a(i,j) a(i,j+1)][oz=0 l=30u+l+150u+50u w=0u h=0u]      
   ]  
]
gapcloser p1 [a(1)][l=100u w=20u oz=0 gap=10u h=4u]

