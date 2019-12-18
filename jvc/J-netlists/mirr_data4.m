%equivalent mass of comb array attachment
h=50e-6; 
N=70; %number fingers and holes
gap=3e-6;
fw=3.2e-6; %finger width
LG=(2*gap+2*fw); %spacing
w2G=(30e-6 - 12e-6)/2; %whoriz
w1G=6e-6; %wvert
wG=30e-6; %wperf
lg=30e-6; %finger length
B=N*lg*fw*h; %vol fingers
A=(2*LG*w2G*h+wG*w1G*h)*N; %vol arm
L=4*140e-6; %length support
w=30e-6; %wperf
w1=8e-6; %wvert
w2=10e-6; %whoriz
C=2*L*h*w2 + N*w1*w*h; %vol support
vol=8*(A+B) + C

LW=vol/h;
L=4*140e-6
W=LW/L
h
h*W*L
