function [sigma,s,totcharge,charge,f,F,sig]=electro3d3(N,gap,volt,a)
k=0;
for i=1:N
   for j=1:N      
      k=k+1;
      n{k}=[0 0 -1]';
      r{k}=[i j gap]';
      da(k)=pi*a^2;
      V(k,1)=volt(1);
   end
end
for i=1:N
   for j=1:N
      k=k+1;
      n{k}=[0 0 1]';
      r{k}=[i j 0]';
      da(k)=pi*a^2;
      V(k,1)=volt(2);
   end
end
fpe=4*pi*8.854e-12;
%V=sum(S*sigma/fpe) where S=da/(norm(rj-ri)).
epsillon=1e-30;S=zeros(2*N,2*N);
for i=1:k
   for j=1:k
      if (i~=j)
         S(i,j)=da(j)/norm(r{i}-r{j}+epsillon,2);
s1=da(j)/norm(r{i}-r{j}+epsillon,2)
s2=2*pi*(sqrt(norm(r{i}-r{j}+epsillon,2)^2+a^2)-norm(r{i}-r{j}+epsillon,2))
s1s2=s1-s2
sigma1=20/s1*fpe
sigma2=20/s2*fpe
V1=sigma1*da(j)/fpe/norm(r{i}-r{j}+epsillon,2)

%S(i,j)=2*pi*(sqrt(1+0.5^2)-1)
      else
         S(i,j)=2*pi*a; %approximation to the square da. Supposed to be circumference!
%'3',         4*sqrt(da(j))
      end
   end
end
sigma=(S\V)*fpe; %surface_ charge calculation
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
%mesh(sig);
