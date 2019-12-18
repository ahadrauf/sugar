function [residual] = elliptic_p_of_Fphi0_minus_FphiB(p,thetaB2,L3k,L1k)
    phi0 = acos(L3k/2/p);
    phiB = asin(sin(thetaB2)/p);
    FphiB = quadl('elliptic',0,phiB,[],[],p,1);
    Fphi0 = quadl('elliptic',0,phi0,[],[],p,1);
    residual = norm(Fphi0 - FphiB - L1k);

    
    