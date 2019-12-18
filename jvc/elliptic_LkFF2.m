function [residual] = elliptic_LkFF2(p,Lk,p1)
    phi1 = asin(sqrt(1/p^2 - p1));
    F1 = quadl('elliptic',0,phi1,[],[],p,1); %incomplete elliptic integral of the first kind
    F2 = quadl('elliptic',0,pi/4,[],[],p,1);

    %residual = real( Lk/p - F1 + F2) + imag( Lk/p - F1 + F2);
    
    residual = norm( Lk/p - F1 + F2);
        
    
    