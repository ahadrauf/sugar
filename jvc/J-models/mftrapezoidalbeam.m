%p.angle=90;p.l2=150e-6;net=cho_load('mftriangle.m',p);figure(1);cho_display(net);
%net=cho_load('mftriangle.m');figure(1);cho_display(net);
uses mumps2.net
param l, wtop, wbottom, h

%straight-a-ways
beam3d p1 [A1t A2t][l=l/4 w=wtop h=H/3] 
beam3d p1 [A1m A2m][l=l/4 w=(wtop+wbottom)/2 h=H/3] 
beam3d p1 [A1b A2b][l=l/4 w=wbottom h=H/3] 

beam3d p1 [A2t A3t][l=l/4 w=wtop h=H/3] 
beam3d p1 [A2m A3m][l=l/4 w=(wtop+wbottom)/2 h=H/3] 
beam3d p1 [A2b A3b][l=l/4 w=wbottom h=H/3] 

beam3d p1 [A3t A4t][l=l/4 w=wtop h=H/3] 
beam3d p1 [A3m A4m][l=l/4 w=(wtop+wbottom)/2 h=H/3] 
beam3d p1 [A3b A4b][l=l/4 w=wbottom h=H/3] 

beam3d p1 [A4t A5t][l=l/4 w=wtop h=H/3] 
beam3d p1 [A4m A5m][l=l/4 w=(wtop+wbottom)/2 h=H/3] 
beam3d p1 [A4b A5b][l=l/4 w=wbottom h=H/3] 

%connectors
beam3d p1 [A1t A1m][l=H/3 w=wtop h=l/4/2 oy=-90] 
beam3d p1 [A1m A1b][l=H/3 w=wbottom h=l/4/2 oy=-90] 

beam3d p1 [A2t A2m][l=H/3 w=wtop h=l/4 oy=-90] 
beam3d p1 [A2m A2b][l=H/3 w=wbottom h=l/4 oy=-90] 

beam3d p1 [A3t A3m][l=H/3 w=wtop h=l/4 oy=-90] 
beam3d p1 [A3m A3b][l=H/3 w=wbottom h=l/4 oy=-90] 

beam3d p1 [A4t A4m][l=H/3 w=wtop h=l/4 oy=-90] 
beam3d p1 [A4m A4b][l=H/3 w=wbottom h=l/4 oy=-90] 

beam3d p1 [A5t A5m][l=H/3 w=wtop h=l/4/2 oy=-90] 
beam3d p1 [A5m A5b][l=H/3 w=wbottom h=l/4/2 oy=-90] 

