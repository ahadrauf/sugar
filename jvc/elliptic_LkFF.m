function [residual] = elliptic_LkFF(p,Lk,p1,p2)
    phi1 = acos(p1/p);
    phi2 = asin(p2/p);
    F1 = quadl('elliptic',0,phi1,[],[],p,1); %incomplete elliptic integral of the first kind
    F2 = quadl('elliptic',0,phi2,[],[],p,1);
    residual = norm( Lk - (F1 - F2) );
        