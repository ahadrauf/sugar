function [x,y,rz,Fx,Fy,Mz]=nl_1(rz,fr);

%Given rz, calculate x, y, fx, fy, m, p, k, phib
    phi0=pi/2 *fr; 
    thetaB=pi/2; 
    psi0=rz; 

%test case for beam shapped like a hanger hook    
    %Mz=7e-9;
    %Fy=-220e-6;
    L3=-3.18181818181818e-5
    k=31622.7766016838
    phiB=1.5707963267949
    phi0=2.799898473665
    psi0=0.340909090909091
    thetaB=pi/2
    p=L3*k/2/cos(phi0)
L3_=2*p*cos(phi0)/k
p
cosphi0=cos(phi0)
k
%integration results    
    [Fphi0,FphiB,Ephi0,EphiB,phiB,p,theta0]=NL_integrations1(phi0,thetaB,psi0);
%    [Fphi0,FphiB,Ephi0,EphiB]=NL_integrations1(phi0,thetaB,psi0);
    
    E=165e9;
    I=(2e-6)^4/12;
    L1=50e-6;
k=(Fphi0-FphiB)/L1; 
    P=k*k*E*I;
    x=2*p*(cos(phiB)-cos(phi0))/k; %lateral deflection at tip
    y=(2*EphiB-FphiB-2*Ephi0+Fphi0)/k; %transverse deflection at tip
    dx=L1-x; %beam shortening
    Fx=P*cos(thetaB);
    Fy=P*sin(thetaB); 
L3=2*p*cos(phi0)/k; %moment arm
    Mz=L3*Fy;
    
%plot beam    
if 1    
    i=0; s=0;
    for phi=phiB:(phi0-phiB)/50:phi0    
        i=i+1;
        X(i)=2*p*(cos(phiB)-cos(phi))/k; %lateral deflection at tip
        [Fphi,FphiB,Ephi,EphiB,phiB,p,theta0]=NL_integrations1(phi,thetaB,psi0);
%        [Fphi,FphiB,Ephi,EphiB]=NL_integrations1(phi,thetaB,psi0);
        Y(i)=(2*EphiB-FphiB-2*Ephi+Fphi)/k; %transverse deflection at tip
        if i>1
            s=s+sqrt((Y(i)-Y(i-1))^2 + (X(i)-X(i-1))^2);
        else
            s=s+sqrt((Y(i))^2 + (X(i))^2);
        end            
    end
    figure(1);
    plot(X,-Y,'g');
end

s
%p.fx=0;p.fy=220e-6;p.m=-7e-9; net=cho_load('can06.m',p); q=cho_dc(net); figure(1); cho_display(net,q); qx=q(lookup_coord(net,'B','x')),qy=q(lookup_coord(net,'B','y')),qrz=q(lookup_coord(net,'B','rz'))

%beam like a hanger hook 
if 0
    Mz=7e-9;
    Fy=-220e-6;
    L3=Mz/Fy
    k=sqrt(abs(Fy/(E*I)))
    phiB=pi/2
    phi0=fzero(@fn_phi0,pi/2*-0.9,[],phiB,p,L1,k)
    %Mz=7e-9;
    %Fy=-220e-6;
    %L3=-3.18181818181818e-5
    %k=31622.7766016838
    %phiB=1.5707963267949
    %phi0=2.799898473665
    %psi0=0.340909090909091
end       

function err=fn_phi0(phi0,phiB,p,L1,k)
    EphiB=quadl(@ellipticintegrand,0,phiB,[],[],p,2); %elliptic integral of the second kind 
    FphiB=quadl(@ellipticintegrand,0,phiB,[],[],p,1); %elliptic integral of the fist kind 
    Ephi0=quadl(@ellipticintegrand,0,phi0,[],[],p,2); %elliptic integral of the second kind 
    Fphi0=quadl(@ellipticintegrand,0,phi0,[],[],p,1); %elliptic integral of the first kind 
    err=k-(Fphi0-FphiB)/L1;   
    
function F=ellipticintegrand(x,p,kind)
%%By: jvclark - spring 2002
    if kind==1 %elliptic integral of the first kind
        F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
    else %elliptic integral of the second kind
        F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
    end

