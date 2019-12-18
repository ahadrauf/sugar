subnet combarray [a(1) b][N=* w=* h=* l=* angle=*] 
[
    wb = 2 * l * tan(angle)
    for j=1:N-1
    [
        W = w + (wb-w)*(j)/(N)
        beam3d p1 [a(j) a(j+1)][l=l/N w=W h=h ] 
    ]
    j = N
    W = w + (wb-w)*(j)/(N)
    beam3d p1 [a(j) b][l=l/N w=W h=h ] 
]

   
 