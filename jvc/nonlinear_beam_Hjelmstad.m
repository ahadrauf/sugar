function nonlinear_beam_Hjelmstad

%example
EI=1e6;     %bending modulus
GA=1e6;     %shear modulus
EA=1e3;     %bulk modulus
L=10;       %length
lambda=[0,1,2,3,4]*pi/2; %Moment factor
Mo=6.28;    %EI/L;    %moment
radius_of_curvature=Mo/EI;
Vo=0;       %lateral force
Ho=0;       %axial force
tol=0.1e-7; %tolerance
a=10;       %step length, arc length parameter

%problem parameters
tol=tol;    %convergence tolerance, solution tolerance
alpha=a;    %step size, arc length parameter
xlength=10; %column length
maxsteps=10/5;%number of steps, number of load steps
maxit=40;   %max number of iterations
inter=10;   %number of integration intervals
npts=20;    %number of output points, , number of points in plotted shape
nnstep=10;  %output interval, output shape every N points, N
nbasis=6;   %number of basis functions
ndm=3*nbasis+1;
yes=1;no=0;

D(1)=EI;    %bending modulus
D(2)=GA;    %shear modulus
D(3)=EA;    %bulk modulus
Cm=Mo;      %concentrated end loads Ml - moment
Cq=Vo;      %concentrated end loads Vl - lateral force
Cp=Ho;      %concentrated end loads Hl - axial force
dmo(1)=0;   %distributed moment loads m1
dmo(2)=0;   %distributed moment loads m2
dmo(3)=0;   %distributed moment loads m3
dqo(1)=0;   %distributed transverse loads q1
dqo(2)=0;   %distributed transverse loads q2
dqo(3)=0;   %distributed transverse loads q3
dpo(1)=0;   %distributed axial loads p1
dpo(2)=0;   %distributed axial loads p2
dpo(3)=0;   %distributed axial loads p3

A=zeros(ndm,ndm);
b=zeros(ndm,1);
x=zeros(ndm,1);
xo=zeros(ndm,1);

%initialize integration increment, change intervals to Simpson points
dz = 1/(2*inter);
npoints = 2*inter + 1;

%initialize values for load step zero, set first load level
x(ndm) = alpha;

fprintf('\n\nFully nonlinear beam analysis:');
fprintf('\n Convergence tolerence =         %0.4g',tol);
fprintf('\n arc length parameter =          %0.4g',alpha);
fprintf('\n number of load steps =          %0.4g',maxsteps);
fprintf('\n maximum number of iterations =  %0.4g',maxit);
fprintf('\n no. of integration intervals =  %0.4g',inter);
fprintf('\n no. of points in plotted shape= %0.4g',npts);
fprintf('\n output shape every N points, N= %0.4g',nnstep);
fprintf('\n number of basis functions =     %0.4g',nbasis);

fprintf('\n\nBeam properties, end load and dist. load amplitudes');
fprintf('\n column length: \t%0.4g',xlength);
fprintf('\n EA GA EI     : \t%0.4g\t\t%0.4g\t\t%0.4g',D(1),D(2),D(3));
fprintf('\n Mo Qo Po     : \t%0.4g\t\t%0.4g\t\t\t%0.4g',Cm,Cq,Cp);
fprintf('\n mo1 qo1 po1  : \t%0.4g\t\t\t%0.4g\t\t\t%0.4g',dmo(1),dqo(1),dpo(1));
fprintf('\n mo2 qo2 po2  : \t%0.4g\t\t\t%0.4g\t\t\t%0.4g',dmo(2),dqo(2),dpo(2));
fprintf('\n mo3 qo3 po3  : \t%0.4g\t\t\t%0.4g\t\t\t%0.4g',dmo(3),dqo(3),dpo(3));

fprintf('\n\n Step \t Load \t\t u(L) \t\t\t w(L) \t\t theta(L) \t\t NU \t\t ||b||');

%compute maxsteps points along the equilibrium path
%2
for nstep=1:maxsteps 
    
    %perform Newton interation at each load step
    nu = 0;
    goto1 = yes;
    while (goto1)%goto1
        %one
        nu = nu+1;    
    
        %execute numberical integration of Hessian and residual components
        b=zeros(ndm,1);
        A=zeros(ndm,ndm);        
        z=0;
        for mpoint=1:npoints
            [factor]=simpson(mpoint,npoints,dz,xlength);
            [A,b]=fcn(D,x,z,xlength,factor,ndm,nbasis,dmo,dqo,dpo,A,b);
            z=z+dz;
        end
        
        %add end load terms to the residual and coefficient matrix 
        for i=1:nbasis
            [h,dh]=basis(i,1,xlength);
            mm=3*(i-1);
            b(mm+1)=b(mm+1) - h*Cp*x(ndm);
            b(mm+2)=b(mm+2) - h*Cq*x(ndm);
            b(mm+3)=b(mm+3) - h*Cm*x(ndm);
            A(mm+1,ndm)=A(mm+1,ndm) - h*Cp;
            A(mm+2,ndm)=A(mm+2,ndm) - h*Cq;
            A(mm+3,ndm)=A(mm+3,ndm) - h*Cm;
        end
        
        %add arc-length constraint terms to Hessian and residual
        for k=1:ndm
            b(ndm)=b(ndm) + (x(k) - xo(k))^2;
            A(ndm,k)=2*(x(k) - xo(k));
        end
        b(ndm)=b(ndm) - alpha^2;
        
        %compute norm of residual for convergence test
        test = 0;
        for k=1:ndm
            test=test+b(k)^2;
        end
        test = sqrt(test);
        
        %update state vector        
%        A=inv(A);
        [A]=invert(A,ndm,ndm);
        for j=1:ndm
            for k=1:ndm
                x(j) = x(j) - A(j,k)*b(k);
            end
        end
        
        %test for convergence, is successful output values
        if ((test>tol) & (nu<maxit))
            goto1=yes;
        else
            goto1=no;
        end
    end%goto1
    
    if (nu>=maxit)
        'Iteration limit exceeded'
        break;            
    end
    results(x,xlength,nstep,nu,test,ndm,nbasis,npts,nnstep);
        
    %update values for previous converged state and guess at next state
    for j=1:ndm
        temp=xo(j);
        xo(j)=x(j);
        x(j)=2*x(j) - temp;
    end
end %2

'Maximum number of steps exhausted'

%======================================
%--------------------------------------

%compute contribution to A and b at current integration point
function [A,b]=fcn(D,x,z,xlength,factor,ndm,nbasis, dmo,dqo,dpo,A,b)
GB=zeros(4,4);G=zeros(4,3);BGB=zeros(3,3);
%compute displacements and derivative
du=0;dw=0;dtheta=0;theta=0;
for i=1:nbasis
    [h,dh]=basis(i,z,xlength);
    du=du+x(3*i-2)*dh;
    dw=dw+x(3*i-1)*dh;
    dtheta=dtheta+x(3*i)*dh;
    theta=theta+x(3*i)*h;
end
ct=cos(theta);
st=sin(theta);

%compute axial strain, shear strain, and curvature
epsilon = dw*st + (1 + du)*ct - 1;
beta = dw*ct - (1 + du)*st;
curv = dtheta;

%compute axial force, shear force, bending moment, and other forces
bend = D(1)*curv;
shear = D(2)*beta;
axial = D(3)*epsilon;
Hor=axial*ct - shear*st;
Ver=axial*st + shear*ct;
Xi = (1 + du)*Hor + dw*Ver;
Yi = dw*Hor - (1 + du)*Ver;

%compute components of E'DE+G store in matrix G
G(1,1) = D(2)*st*st + D(3)*ct*ct;
G(1,2) = ct*st*(D(3)-D(2));
G(1,3) = D(2)*st*(1+epsilon) + D(3)*ct*beta - Ver;
G(1,4) = 0;
G(2,2) = D(2)*ct*ct + D(3)*st*st;
G(2,3) = D(3)*st*beta - D(2)*ct*(1 + epsilon) + Hor;
G(2,4) = 0;
G(3,3) = D(2)*(1 + epsilon)^2 + D(3)*beta^2 - Xi;
G(3,4) = 0;
G(4,4) = D(1);

%compute the rest of E'DE+G by symmetry
for i=1:3
    for j=i+1:4
        G(j,i)=G(i,j);
    end
end

%Form stiffness matrix K and store it in matrix A
for i=1:nbasis
    [hi,dhi]=basis(i,z,xlength);
    for j=1:nbasis
        [hj,dhj]=basis(j,z,xlength);

        %compute B'(E'DE+G)B noting the sparse structure of B
        for k=1:4
            GB(k,1)=dhj*G(k,1);
            GB(k,2)=dhj*G(k,2);
            GB(k,3)= hj*G(k,3) + dhj*G(k,4);
        end
        for k=1:3
            BGB(1,k)=dhi*GB(1,k);
            BGB(2,k)=dhi*GB(2,k);
            BGB(3,k)= hi*GB(3,k) + dhi*GB(4,k);
        end
        
        %assemble the result into the matrix A
        for m=1:3
            mm=3*(i-1);
            for n=1:3
                nn=3*(j-1);
                A(mm+m,nn+n)=A(mm+m,nn+n) + BGB(m,n)*factor;
            end
        end
    end

    %Form integral part of resisdual force and assemble into matrix b
    mm=3*(i-1);
    [dm,dq,dp]=applied(z,x(ndm),1,dmo,dqo,dpo);
    b(mm+1) = b(mm+1) + (dhi*Hor - hi*dp)*factor;
    b(mm+2) = b(mm+2) + (dhi*Ver - hi*dq)*factor;
    b(mm+3) = b(mm+3) + (dhi*bend + hi*(Yi - dm))*factor;

    %Form the integral part of load factor part of matrix A
    [dm,dq,dp]=applied(z,x(ndm),2,dmo,dqo,dpo);
    A(mm+1,ndm) = A(mm+1,ndm) - hi*dp*factor;
    A(mm+2,ndm) = A(mm+2,ndm) - hi*dq*factor;
    A(mm+3,ndm) = A(mm+3,ndm) - hi*dm*factor;
end

%-----------------------

%evaluate ith basis function h and derivative dh
function [h,dh]=basis(i,z,xlength)
n = mod(i,2);
m = i/2 + n;
a = (2*m-1)*2*atan(1); 
if (n==0)
    h = sin(a*z);
    dh = a*cos(a*z)/xlength;
else
    h = 1 - cos(a*z);
    dh = a*sin(a*z)/xlength;
end

%------------------------

%evaluate the distributed load functions at point z
function [dm,dq,dp]=applied(z,CLf,n,dmo,dqo,dpo)

%compute the nominal vaules of the applied forces at point z
f = 1 - z;
if n==1
    dm = dmo(1) + (dmo(2) + dmo(3)*f)*CLf;
    dq = dqo(1) + (dqo(2) + dqo(3)*f)*CLf;
    dp = dpo(1) + (dpo(2) + dpo(3)*f)*CLf;
end
if n==2
    dm = dmo(2) + dmo(3)*f;
    dq = dqo(2) + dqo(3)*f;
    dp = dpo(2) + dpo(3)*f;
end

%-----------------------------------

%print results of current step to various files
function results(x,xlength,nstep,nu,test,ndm,nbasis,npts,nnstep)
yes=1;no=0;
if (mod(nstep,nnstep)==0)
    PLOT=yes;
    nplot=npts;
else
    PLOT=no;
    nplot=1;
end

%write coefficients for current step
if PLOT 
%    fprintf('\nLoad Step = %0.4g \t\t Load factor = %0.4g \t\t Iterations = %0.4g \t\t Norm of Residual = %0.4g',nstep,x(ndm),nu,test);
    for i = 1:nbasis
        ii=3*(i-1);
%        fprintf('\n i = %0.4g \t\t a(i) = %0.4g \t\t b(i) = %0.4g \t\t c(i) = %0.4g',i,x(ii+1),x(ii+2),x(ii+3));
   end
end

%compute and print current geometry of beam
if PLOT
%    fprintf('\n Step = %0.4g \t\t Load = %0.4g',nstep,x(ndm));
end
z=0;
dz = 1/nplot;
for ii = 1:nplot+1
    for j=1:3
        disp(j) = 0;
        for i=1:nbasis
            [h,dh]=basis(i,z,xlength);
            disp(j) = disp(j) + h*x(3*(i-1)+j);
        end
    end
    if PLOT
%        fprintf('\n z*xlength=%0.4g,  z*xlength+disp(1)=%0.4g,  disp(2)=%0.4g',z*xlength,z*xlength+disp(1),disp(2));
    end
    z=z+dz;
end

fprintf('\n   \t%d \t\t%0.4g \t\t%0.4g \t%0.4g \t%0.4g \t\t%0.4g \t%0.4g',nstep,x(ndm),disp(1),disp(2),disp(3),nu,test);

%-------------------------------

%evaluate Simposn intrgration weight factors
function [factor]=simpson(m,npoints,dz,xlength)
c = xlength*dz/3;
n = mod(m,2);
if ((m==1) | (m==npoints))
    factor = c;
elseif n==0
    factor = 4*c;
else
    factor = 2*c;
end

%----------------------------
 
%invert matrix A (nmax,nmax) when array dimension is ndm
function [a]=invert(a,nmax,ndm)
for n=1:nmax
    d = a(n,n);
    for j = 1:nmax
        a(n,j) = -a(n,j)/d;
    end
    for i = 1:nmax
        if n~=i
            for j = 1:nmax
                if n~=j
                    a(i,j)=a(i,j) + a(i,n)*a(n,j);
                end
            end
            a(i,n) = a(i,n)/d;
        end
        a(n,n) = 1/d;
    end
end

