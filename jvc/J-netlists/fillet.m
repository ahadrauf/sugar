%Fillet
%fillet p1 [a(1) b] [l=100u N=10 angle=10*pi/180 w=2u h=H ] 

subnet fillet [a(1) b][N=* w=* h=* l=* wb=*] 
[
    for j=1:N-1
    [
        W = w + (wb-w)*(j)/(N)
        beam_link_grav p1 [a(j) a(j+1)][l=l/N w=W h=h ] 
    ]
    j = N
    W = w + (wb-w)*(j)/(N)
    beam_link_grav p1 [a(j) b][l=l/N w=W h=h ] 
]

   
 