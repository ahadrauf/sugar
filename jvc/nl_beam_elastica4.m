function nl_beam_elastica4

W = 100e-6; %beam width. y axis.
L = 100e-6; %beam length. x axis.
H = 2e-6; %beam layer thickness. z axis.
E = 165e9; %Young's modulus.
I = W*H^3 / 12; %moment of inertial about y.
EI = E*I; %flexural rigidity
F = 10e-6; %tip force
M = 10e-8; %tip moment
e = M/F; %tip moment arm
k = F/E/I; %the k parameter. units = 1/L.
p_guess = 0.6; %0 < p < 1
p1 = e*k/2;
p2 = 1/sqrt(2);



clear all
fprintf('For e2 ');
%find k and e2
F = 1
EI = 10000
k = sqrt(F/EI)
M = 50*2
M2P_lt_2EI = M*M/F < 2*EI
e2 = M/F 
%find p = p(e,k)
L = 100
Lk = L*k;
p1 = (e2*k/2)^2
%p = fminbnd('elliptic_LkKFF', 0,2,optimset('tolx',eps), L,k, p1) 
p = fminsearch('elliptic_LkKFF', 0.7 ,optimset('tolx',eps), L,k, p1) 
%other params    
h = 2/p/k 
ymin = h*sqrt(1-p*p)
phi1 = asin(sqrt(1/p^2 - p1));
F1 = quadl('elliptic',0,phi1,[],[],p,1); %incomplete elliptic integral of the first kind
F2 = quadl('elliptic',0,pi/4,[],[],p,1);
K  = quadl('elliptic',0,pi/2,[],[],p,1); %complete elliptic integral of the first kind
S = p/k*( 2*K - F1 - F2)
thetaD = phi1*2