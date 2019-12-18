function [sigma,s,totcharge,charge,f,F,sig]=electro3d2(N,gap,volt)

k=0;
for i=1:N
   for j=1:N      
      k=k+1;
      n{k}=[0 0 -1]';
      r{k}=[i j gap/2]';
      da(k)=pi*0.5*0.5;
      V(k,1)=volt(1);
   end
end
for i=1:N
   for j=1:N
      k=k+1;
      n{k}=[0 0 1]';
      r{k}=[i j -gap/2]';
      da(k)=pi*0.5*0.5;
      V(k,1)=volt(2);
   end
end

%V=sum(S*sigma/(4*pi*e0)) where S=da/(abs(rj-ri)).
epsillon=1e-30;S=zeros(2*N,2*N);ii=0;jj=0;
for i=1:k
   for j=1:k
      if (i~=j)
         S(i,j)=da(j)/norm(r{i}-r{j}+epsillon,2);
      else
         S(i,j)=4*sqrt(da(j)); %approximation to a square da. Supposed to be 4*circumference!
      end
   end
end
sigma=(S\V)*4*pi*8.854e-12; %surface_ charge calculation
s=sum(sigma);

%forces per patch
f=zeros(3,k);
for i=1:k
   f(1:3,i)=sigma(i)^2*da(i)*n{i}/2/8.854e-12;
end

%forces per surface
%F=zeros(3,2);
F=zeros(3,1);
for i=1:k/2
%   F(1:3,1+(i>k/2))=F(1:3,1+(i>k/2))+f(1:3,i);
   F(1:3,1)=F(1:3,1)+f(1:3,i);
end
   
%charge per surface
charge=zeros(k,1);
for i=1:k
   charge(i)=sigma(i)*da(i);
end
totcharge=sum(charge(1:k/2));

%capacitance of surface i due to j
%C=zeros(k,k); 
%for i=1:k
%   for j=1:k
%      C(i,j)=abs(charge(i)/(V(i)-V(j)+epsillon)); %capacitance      
%   end
%end
I=0;
figure(1);hold on;
for i=1:N
   for j=1:N
      I=I+1;
      %      plot3(i,j,sigma(I));
      sig(i,j)=sigma(I);
   end
end
mesh(sig);
