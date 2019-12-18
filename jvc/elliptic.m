%function F=elliptic(x,p,kind)
% *Complete* elliptic integral of the first kind is K(pi/2,p)
% *Complete* elliptic integral of the first kind is E(pi/2,p)
%The Matlab function for the complete elliptic integrals are [K,E] = ellipke(p*p), where the input modulus is squared.
%The equivalent 2nd kind is [E]=quadl('elliptic',0,pi/2,[],[],p,2).

function F=elliptic(x,p,kind)
if kind==1
   F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind   
else
   F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
end

