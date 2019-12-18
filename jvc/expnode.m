function [spx,spy]=expnode(X,Y,Z)
Z=real(Z); 
[gx,gy]=ginput(4);
[lenx,leny]=size(Z);

%determine x slope
y=(gy(1)+gy(2))/2; %average out y
ii=0;

i1=ceil(gx(1)/max(max(X))*lenx);
i2=floor(gx(2)/max(max(X))*lenx);   
j1=ceil(gy(1)/max(max(Y))*leny);
j2=floor(gy(2)/max(max(Y))*leny);   
yave=(Y(i1,j1)+Y(i1,j1))/2;
resolution=max(abs(i2-i1),abs(j2-j1));
inci=abs(i2-i1)/resolution;
incj=abs(j2-j1)/resolution;
i=i1; j=j1;
for k=1:resolution   
   i=round(i+(k-1)*inci);
   j=round(j+(k-1)*incj);
   z(k)=Z(i,j);
   x(k)=abs(sqrt(X(i,j)^2+Y(i,j)^2)-sqrt(X(i1,j1)^2+Y(i1,j1)^2));
end   
   
   if Z(
x=X(

[lenx,leny]=size(Z);
%along x
ii=0;
xi=ceil(lenx/2);
yi=ceil(leny/2);
xi=60;
for i = 1:leny
%for i = 40:60
   if Z(xi,i)~=0
      ii =ii+1;
      spx(ii)=X(xi,i);
      spy(ii)=Z(xi,i);
   end
end

slopealongx = mmspder(spx,spy,Y(xi,yi))

%=========================
function z=mmspder(x,y,xi)

if nargin==3
   pp=spline(x,y);
else
   pp=x;
end
[br,co,npy,nco]=unmkpp(pp); %take apart pp
sf=nco-1:-1:1;
dco=sf(ones(npy,1),:).*co(:,1:nco-1);
ppd=mkpp(br,dco);
if nargin==1
   z=ppd;
elseif nargin==2
   z=ppval(ppd,y);
else
   z=ppval(ppd,xi);
end

