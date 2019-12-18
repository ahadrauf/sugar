function wm
format compact

k=4;
m1=10;
w1=sqrt(k./m1);
m2=20;
w2=sqrt(k./m2);
m3=30;
w3=sqrt(k./m3);
m4=40;
w4=sqrt(k./m4);


w12=w1-w2;
w23=w2-w3;

for k=1:5
    k
    DM43 = abs(k/w4^2 - k/w3^2) 
    DM12 = abs(k/w1^2 - k/w2^2)
end
