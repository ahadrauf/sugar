%function [Fphi0,Fphib,Ephi0,Ephib,phib,p,theta0]=NL_integrations1(phi0,thetaB,psi0) 

function [Fphi0,FphiB,Ephi0,EphiB,phiB,p,theta0]=NL_integrations1(phi0,thetaB,psi0) 
%By: jvclark - spring 2002
    
%parameters
    theta0=thetaB+psi0; %thetaB => angle of force, psi0 => tip angle rz.
    p=sin(theta0/2); %elliptic integral modulus
    phiB=asin(sin(thetaB/2)/p); %angle associated with beam length and moment
    
%elliptic integration
    EphiB=quadl(@ellipticintegrand,0,phiB,[],[],p,2); %elliptic integral of the second kind 
    FphiB=quadl(@ellipticintegrand,0,phiB,[],[],p,1); %elliptic integral of the fist kind 
    Ephi0=quadl(@ellipticintegrand,0,phi0,[],[],p,2); %elliptic integral of the second kind 
    Fphi0=quadl(@ellipticintegrand,0,phi0,[],[],p,1); %elliptic integral of the first kind 
        
function F=ellipticintegrand(x,p,kind)
%%By: jvclark - spring 2002
    if kind==1 %elliptic integral of the first kind
        F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
    else %elliptic integral of the second kind
        F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
    end

