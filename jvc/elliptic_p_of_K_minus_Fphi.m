function [residual] = elliptic_p_of_K_minus_Fphi(p,L1k)
    [K, E] = ellipke(p*p);
    phiB = asin(1/p/sqrt(2));
    FphiB = quadl('elliptic',0,phiB,[],[],p,1);
    residual = norm( K - FphiB - L1k);

    