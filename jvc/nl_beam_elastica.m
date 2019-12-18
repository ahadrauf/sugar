%nl_beam_elastica(length, width, thickness, Young's modulus, force, moment)

function nl_beam_elastica(L,W,H,E,F,M)
test = 1;
if nargin==0
    W = 2e-6; %beam width. y axis.
    L = 100e-6; %beam length. x axis.
    H = 2e-6; %beam layer thickness. z axis.
    E = 165e9; %Young's modulus.
    F = 1e-6; %tip force
    M = 0; %tip moment
end

I = W*H^3 / 12; %moment of inertial about y.
EI = E*I; %flexural rigidity

if F == 0
    F = M;
    e = 1;
else
    e = M/F;
end

%if a test is to be made, EI=10000, F=1, M=50, & L=100, produces p=0.974.
if test
    fprintf('This is a test: EI=10000, F=1, M=50, & L=100, should produce p=0.974.\n');
    EI=10000; F=1; M=50; e=M/F; L=100; k=sqrt(F/EI); p1=e*k/2; %=> p=0.974
end

%determin p parameter 
%p_guess = 0.6; %0 < p < 1
k = F/EI; %the k parameter. units = 1/L.
p1 = e*k/2;
p2 = 1/sqrt(2);
%p = fzero('elliptic_LkFF', p_guess, [], L*k, p1, p2); 
p = fminbnd('elliptic_LkFF', eps,1-eps, [], L*k, p1, p2); 

%check for undulating elastica or nodal elastica
if (M^2/F < 2*EI) & (p < 1) & (p > 0)
    undulating_elastica = 1;
else 
    undulating_elastica = 0;
    fprintf('Warning: The force and moment chosen produces nodal elastica => tip > pi/2.\n');
end

if undulating_elastica
    %0 < tip < pi/2.
    fprintf('This deflection type is undulating elastica. \n');    
    phiC = acos(e*k/2/p);
    thetaC = 2 * asin(p*sin(phiC));
    phiB = asin(1/p/sqrt(2));
    phi = phiC;
    ii = 0;
    for phi = phiB : (phiC-phiB)/100 : phiC
        ii = ii + 1;
        EphiB = quadl('elliptic',0,phiB,[],[],p,2);
        Ephi = quadl('elliptic',0,phi,[],[],p,2);
        FphiB = quadl('elliptic',0,phiB,[],[],p,1);
        Fphi = quadl('elliptic',0,phi,[],[],p,1);
        x = 2*p*(cos(phiB) - cos(phi))/k;
        y = (2*EphiB - FphiB + Fphi - 2*Ephi)/k;

        X(ii) = x;
        Y(ii) = y;
    end
    figure(1); plot(X,Y); grid on; axis equal;
end

if test
    p
end
Y_nonlin = max(Y)
Y_linear = (F*L^3)/(3*EI)
L 
F
EI
