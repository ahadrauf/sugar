%tic; b_=10e-9;c_=100;A_=2e-6;V0_=15; [V,x,y]=surfaceroughness_V2d_4(b_,c_,A_,V0_);toc

function [C]=surfaceroughness_V2d_4

c_ = 10;
A_ = 1e-6;
V0_ = 15;
Xb = 0.5e-6;

%General V(x,y) symbolically
syms A V0 c b x y n;        %symbolic vars
number_summation_terms = 5; %Number of terms in the summation
%V = symsum(2*sinh(n*(A+x))*V0*sin(n*y)*(1+(-1)^(n+1))/(n*pi*sinh(2*n*A)), n, 1, number_summation_terms) + sinh(c*(A+x))*Vb*sin(c*y)/sinh(2*c*A);
Vb = sinh(2*c*A) * V0 * (1 - symsum(2*sinh(n*(A+x))*sin(n*y)*(1+(-1)^(n+1))/(n*pi*sinh(2*n*A)), n, 1, number_summation_terms) )/(sinh(c*(A+x))*sin(c*y));

%b = b_;   %perturbation in V0
c = c_;   %number of cycles / 2
V0 = V0_; %potential contour
A = A_;   %gap / 2

%find Vb, and V(x,y) given physical sinusoidal surface. This assumes that the sinusoidal surface is proportional to a flate sinusoidal potential.
x = A-Xb;           %select particular x
y = pi/2;           %select particular y
Vb_eval = eval(Vb); %numerically evaluate Vb

%potential and electric field anywhere inside the region
syms x y;
V = symsum(2*sinh(n*(A+x))*V0*sin(n*y)*(1+(-1)^(n+1))/(n*pi*sinh(2*n*A)), n, 1, number_summation_terms) + sinh(c*(A+x))*Vb_eval*sin(c*y)/sinh(2*c*A);
Ex = -diff(eval(V),'x');    %x component of E.
Ey = -diff(eval(V),'y');    %y component of E.
E = sqrt(Ex^2 + Ey^2);      %Magnitude of the electric field at (x,y).

lambda_min = 70e-9;         %symmetry only requires a fourth
lambda_max = 200e-9;    
beta_min = 8e-9;            %beta is the amplitude
beta_max = 116e-9;
m=0;        %counter
Q = 0;      %initialize net charge
h = 20e-6;  %layer thickness
L = pi;     %length of the plate in meters
ym = L/2;   %middle y coordinate
permittivity = 8.854187817e-12; %permittivity of free space
resolution = 4;
number_of_segments = 10;
n1 = 0;
for lambda = lambda_min : (lambda_max-lambda_min)/resolution : lambda_max
    for beta = beta_min : (beta_max-beta_min)/resolution : beta_max
        
        %geometric parameters
        c = 2*L/lambda;     %number of cycles / 2
        Xb = beta/2;        %amplitude perturbation in x(y)
        V0 = V0_;           %nominal potential
        A = A_;             %gap / 2
        
        %find Vb, and V(x,y) given physical sinusoidal surface. This assumes that the sinusoidal surface is proportional to a flate sinusoidal potential.
        x = A-Xb;           %select particular x
        y = pi/2;           %select particular y
        Vb_eval = eval(Vb); %numerically evaluate the Vb required to produce that particular x at y 

        %potential and electric field anywhere inside the defined region defined by above geometry
        syms x y;
        V = symsum(2*sinh(n*(A+x))*V0*sin(n*y)*(1+(-1)^(n+1))/(n*pi*sinh(2*n*A)), n, 1, number_summation_terms) + sinh(c*(A+x))*Vb_eval*sin(c*y)/sinh(2*c*A);
        Ex = -diff(eval(V),'x');    %x component of E.
        Ey = -diff(eval(V),'y');    %y component of E.
        x = A-Xb;           %select particular x
        y = pi/2;           %select particular y
        E = sqrt(eval(Ex)^2 + eval(Ey)^2);      %Magnitude of the electric field at (x,y).

        y1 = pi/2;
        Q = 0;
        
        pi/2 + (pi/2+lambda/2)/6 

        abs((pi/2 + (pi/2+lambda/2)/6) - (pi/2+lambda/2))/6 

        (pi/2+lambda/2)
        
        for y2 = pi/2 + (pi/2+lambda/2)/6 : abs((pi/2 + (pi/2+lambda/2)/6) - (pi/2+lambda/2))/6 : (pi/2+lambda/2)
            x1 = A + Xb * sin(c * y1);
            x2 = A + Xb * sin(c * y2);
            s = sqrt( (y2-y1)^2 + (x2-x1)^2 );
            x = (x2 + x1)/2;
            y = (y2 + y1)/2;
            Qn = h * s * permittivity * E
            Q = Q + Qn;
            xn = x; 
            yn = y;
        end
        n1 = n1 + 1;
        C{n1} = {lambda, beta, Xb, c, A, V0, xn, yn, Qn}
    end
end
