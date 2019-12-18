% net=cho_load('nl_trapezoidal.m');[q]=cho_dc(net);cho_display(net,q);

uses mumpsx.net
anchor   p1 [A1m]    [l=10u w=10u h=10u oz=180]
trapezoidalbeam p1 [A1m A5m][l=40u wtop=2u wbottom=6u h=4u oz=0]
%f3d      *  [A5m]    [F=30u oz=90]
