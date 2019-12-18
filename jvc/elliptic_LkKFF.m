function [residual] = elliptic_LkKFF(p,L,k,p1)
    phi1 = pi-asin(sqrt(1/p^2 - p1));
    F1 = quadl('elliptic',0,phi1,[],[],p,1); %incomplete elliptic integral of the first kind
    F2 = quadl('elliptic',0,pi/4,[],[],p,1);
    K  = quadl('elliptic',0,pi/2,[],[],p,1); %complete elliptic integral of the first kind
    n=1;
    S = p/k*((2*n+2)*K - F1 - F2);
    residual = norm( L - real(S));

    %if residual = norm( (100-real(S)) + 0*abs(imag(S)));
        
%    if (abs(100-real(S)) < 1)
%        residual = abs(imag(S));
%    else
%       residual = abs(100-real(S));
%   end
        
 
    
    
   S
