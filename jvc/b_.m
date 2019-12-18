function S=b(x,y,a)
y1=y;
y1=y+1e-30;
a1=(x-a/2)*log(sqrt((x-a/2)^2+y^2));
a2=-(x+a/2)*log(sqrt((x+a/2)^2+y^2));
a3=y*atan((x-a/2)/y1);
a4=-y*atan((x+a/2)/y1)+a;
S=(a1+a2+a3+a4)/2/pi/8.854e-12;

