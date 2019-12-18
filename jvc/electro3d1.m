function [sigma,n]=electro3d1(surfaces,voltages)
%N surfaces => { {X1,Y1,Z1};{X2,Y2,Z2};{X3,Y3,Z3}; ... {XN,YN,ZN}} => a set of surfaces.
%where X1 => [ [x11 x12 x13 x14 ... x1n x11];[x21 x22 x23 x24 ... x2n x21];...;[xm1 xm2 xm3 xm4 ... xmn xm11]]
%where X1,Y1,Z1 are in a format suitable for surfl(X1,Y1,Z1). Each crosssection is a column vector in the matrix X,Y,Z.
%crosssections are assumed to progress clockwise if about positive x
%voltages => [V1 V2 ... VN]

%gather patch areas, normals, voltages, and coordinates. 
X=1;Y=2;Z=3;I=0;
[sectionnodes,numcrosssections]=size(surfaces{surface_}{X}); 
%edge da's. The side edges have four sides and they're planar.
for j=1:numcrosssections-1
   for i=1:sectionnodes-1
      if i==1, u=sectionnodes-1; else u=i-1; end %sectionnodes rap-around effect
      a=[surfaces{surface_}{X}(i,j+1)-surfaces{surface_}{X}(i,j);surfaces{surface_}{Y}(i,j+1)-surfaces{surface_}{Y}(i,j);surfaces{surface_}{Z}(i,j+1)-surfaces{surface_}{Z}(i,j)]; %1st edge vector of da. da has 4 sides. 
      b=[surfaces{surface_}{X}(u,j)-surfaces{surface_}{X}(i,j);surfaces{surface_}{Y}(u,j)-surfaces{surface_}{Y}(i,j);surfaces{surface_}{Z}(u,j)-surfaces{surface_}{Z}(i,j)]; %2nd edge vector of da.
      d=[surfaces{surface_}{X}(u,j+1)-surfaces{surface_}{X}(i,j);surfaces{surface_}{Y}(u,j+1)-surfaces{surface_}{Y}(i,j);surfaces{surface_}{Z}(u,j+1)-surfaces{surface_}{Z}(i,j)]; %diagonal vector of da.
      c=cross(a,b); %normal vector to da.         
      I=I+1;index1(surface_,i,j)=I;index2{I}=[surface_,i,j]; %indices         
      n(1:3,I)=norm(c,2)\c; %normal unit vector of da
      da(I)=2\norm(cross(a,d),2)+2\norm(cross(d,b),2); %area of da
      r(1:3,I)=2\(a+b)+[surfaces{surface_}{X}(i,j);surfaces{surface_}{Y}(i,j);surfaces{surface_}{Z}(i,j)]; %coordinates of da
      V(I,1)=voltages(surface_);
   end %i 
end %j    
%end da's. Assuming end da's are planar polygons.
numtriangles=sectionnodes-3; 
for j=[1,numcrosssections] %end caps
   edgearea=0; edgecoord=0; 
   for k=1:numtriangles
      a=[surfaces{surface_}{X}(k+1,j)-surfaces{surface_}{X}(1,j);surfaces{surface_}{Y}(k+1,j)-surfaces{surface_}{Y}(1,j);surfaces{surface_}{Z}(k+1,j)-surfaces{surface_}{Z}(1,j)]; 
      b=[surfaces{surface_}{X}(k+2,j)-surfaces{surface_}{X}(k+1,j);surfaces{surface_}{Y}(k+2,j)-surfaces{surface_}{Y}(k+1,j);surfaces{surface_}{Z}(k+2,j)-surfaces{surface_}{Z}(k+1,j)]; 
      if j==1, c=cross(b,a); else c=cross(a,b); end %normal vector
      edgearea=edgearea+2\norm(c,2); %sum up triangle areas of the edge         
   end
   for k=1:sectionnodes-1
      edgecoord=edgecoord+[surfaces{surface_}{X}(k,j);surfaces{surface_}{Y}(k,j);surfaces{surface_}{Z}(k,j)]; %coordinates nodes
   end
   I=I+1;index1(surface_,i+1+(j==1),1)=I;index2{I}=[surface_,i+1+(j==1),1]; %indices
   n(1:3,I)=norm(c,2)\c; %normal unit vector of edge da.
   da(I)=edgearea; %area of edge
   r(1:3,I)=(sectionnodes-1)\edgecoord; %average coordinate of the edge nodes.
   V(I,1)=voltages(surface_);
end %j

%V=sum(S*sigma/(4*pi*e0)) where S=da/(abs(rj-ri)).
epsillon=1e-30;S=zeros(I,I);
for i=1:I   
   for j=1:I
      if i~=j
         S(i,j)=da(i)/norm(r(1:3,j)-r(1:3,i)+epsillon,2);
      else
         S(i,j)=4*sqrt(da(i)); %approximation to a square da. Supposed to be 4*circumference!
      end
   end
end
sigma=(S\V)*4*pi*8.854e-12; %surface_ charge calculation
s=sum(sigma);
%forces per patch
f=zeros(3,I);
for i=1:I
   for j=1:I
      if i~=j
%         f(1:3,i)=f(1:3,i)+sigma(j)*da(j)*sigma(i)*da(i)/4/pi/norm(r(1:3,j)-r(1:3,i),2)^2*n(1:3,i);
      end
   end
end

%forces per surface
F=zeros(3,numsurfaces);
for i=1:I
%   F(1:3,index2{i}(1))=F(1:3,index2{i}(1))+f(1:3,i);
end
   
%charge per surface
charge=zeros(I,1);
for i=1:I
%   charge(index2{i}(1))=charge(index2{i}(1))+sigma(i)*da(i);
end

%capacitance of surface i due to j
C=zeros(numsurfaces,numsurfaces); 
for i=1:numsurfaces
   for j=1:numsurfaces
%      C(i,j)=abs(charge(i)/(V(i)-V(j))); %capacitance      
   end
end

