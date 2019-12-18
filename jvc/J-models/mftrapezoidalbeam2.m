%p.angle=90;p.l2=150e-6;net=cho_load('mftriangle.m',p);figure(1);cho_display(net);
%net=cho_load('mftriangle.m');figure(1);cho_display(net);
uses mumps2.net

param l, wtop, wbottom, h
L4=l/4
H3=h/3
H4=h/4
H8=h/8
wtopbot2=(wtop+wbottom)/2

%straight-a-ways
%beam3d p1 [A1t A2t][l=L4 w=wtop h=H3] 
beam3d p1 [A1m A2m][l=L4 w=wtopbot2 h=H3] 
%beam3d p1 [A1b A2b][l=L4 w=wbottom h=H3] 

%beam3d p1 [A2t A3t][l=L4 w=wtop h=H3] 
beam3d p1 [A2m A3m][l=L4 w=wtopbot2 h=H3] 
%beam3d p1 [A2b A3b][l=L4 w=wbottom h=H3] 

%beam3d p1 [A3t A4t][l=L4 w=wtop h=H3] 
beam3d p1 [A3m A4m][l=L4 w=wtopbot2 h=H3] 
%beam3d p1 [A3b A4b][l=L4 w=wbottom h=H3] 

%beam3d p1 [A4t A5t][l=L4 w=wtop h=H3] 
beam3d p1 [A4m A5m][l=L4 w=wtopbot2 h=H3] 
%beam3d p1 [A4b A5b][l=L4 w=wbottom h=H3] 

%connectors
%beam3d p1 [A1t A1m][l=H3 w=wtop h=H8 oy=-90] 
%beam3d p1 [A1m A1b][l=H3 w=wbottom h=H8 oy=-90] 

%beam3d p1 [A2t A2m][l=H3 w=wtop h=H4 oy=-90] 
%beam3d p1 [A2m A2b][l=H3 w=wbottom h=H4 oy=-90] 

%beam3d p1 [A3t A3m][l=H3 w=wtop h=H4 oy=-90] 
%beam3d p1 [A3m A3b][l=H3 w=wbottom h=H4 oy=-90] 

%beam3d p1 [A4t A4m][l=H3 w=wtop h=H4 oy=-90] 
%beam3d p1 [A4m A4b][l=H3 w=wbottom h=H4 oy=-90] 

%beam3d p1 [A5t A5m][l=H3 w=wtop h=H8 oy=-90] 
%beam3d p1 [A5m A5b][l=H3 w=wbottom h=H8 oy=-90] 

