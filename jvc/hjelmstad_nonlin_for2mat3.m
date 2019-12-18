%hjelmstad_nonlin_for2mat3 

function NonlinearBeam 

%   Compute equilibrium path  for Fully Nonlinear Beam      
%      x(3i-2) = a(i), x(3i-1) = b(i), x(3i) = c(i)         
%      a(i) = u coeff, b(i) = w coeff, c(i) = theta coeff   
%      x(ndm) = load factor
%      xo(i) = values of x(i) at last converged state       
%      A = Coefficient matrix, b = right hand side          
%      D = Constitutive properties                          
%      Cm,Cq,Cp Concentrated end forces (M,Q,P)             
%      dmo,dqo,dpo Coefficients for distributed forces      

    nbasis=6;  
    ndm=3*nbasis+1; 

%Read problem parameters
if 0
    tol=input('convergence tolerance = ');    %read(5,1000) tol,alpha,column_length
    alpha=input('arc length parameter = ');
    column_length=input('column length = ');
    maxsteps=input('number of load steps = ');  %read(5,1001) maxsteps,maxit,inter
    maxit=input('maximum number of iterations = ');
    inter=input('number of integration intervals = ');
    npts=input('number of points in plotted shape = ');  %read(5,1001) npts,nnstep
    nnstep=input('output shape every N points, N = ');
    D(1)=input('EA = ');  %read(5,1000) D(1),D(2),D(3)
    D(2)=input('GA = ');
    D(2)=input('EI = ');
    Cm=input('ML = M0 = ');  %read(5,1000) Cm,Cq,Cp
    Cq=input('VL = Q0 = ');
    Cp=input('HL = P0 = ');
    dmo=input('distributed load [m01 m02 m03] = dmo(1,2,3) = '); %read(5,1000) (dmo(i),i=1,3)
    dqo=input('distributed load [q01 q02 q03] = dqo(1,2,3) = '); %read(5,1000) (dqo(i),i=1,3)
    dpo=input('distributed load [p01 p02 p03] = dpo(1,2,3) = '); %read(5,1000) (dpo(i),i=1,3)
else
    tol=0.1e-7;
    alpha=10;
    column_length=10;
    maxsteps=10;
    maxit=40;
    inter=10;
    npts=20;
    nnstep=10;
    D(1)=1000;
    D(2)=0.1e7;
    D(3)=0.1e7;
    Cm=6.28;
    Cq=0;
    Cp=0;
    dmo=[0,0,0];
    dqo=[0,0,0];
    dpo=[0,0,0];
end

fprintf(' Step \t Load \t u(L) \t w(L) \t theta(L) \t   *  \t NU \t || b ||'); 
                               
%Initialize 
    dz = 1/(2*inter);         %integration increment
    npoints = 2*inter + 1;    %change intervals to Simpson points
    [x]=zeros(ndm,1);         %load step
    [xo]=zeros(ndm,1);        %load step
    x(ndm) = alpha;           %first load level
          
%Compute MAXSTEPS points along the Equilibrium Path
for nstep=1:maxsteps   
    nu = 0; %Newton iteration at each load step
    goto_one=1;      
    while goto_one==1      
        nu = nu + 1;        %increment interation counter
        %Numerical integration of Hessian and residual components
        [b]=zeros(ndm,1);   %residual
        [A]=zeros(ndm,ndm); %Hessian
        z = 0;              %Simpson coord
        for mpoint=1:npoints   
          [weight_factor] = simpson(mpoint,npoints,dz,column_length); 
          [A,b] = fcn(A,b,D,x,z,column_length,weight_factor,ndm,nbasis,dmo,dqo,dpo);    
          z = z + dz; %z = z + dz
        end %end do
        %Add end load terms to the residual and coefficient matrix
        for i=1:nbasis   %do i=1,nbasis
          [h,dh]=basis(i,1,column_length,nbasis);    %call basis(i,one,column_length,h,dh,nbasis)
          mm = 3*(i-1);
          b(mm+1) = b(mm+1) - h*Cp*x(ndm);
          b(mm+2) = b(mm+2) - h*Cq*x(ndm);
          b(mm+3) = b(mm+3) - h*Cm*x(ndm);
          A(mm+1,ndm) = A(mm+1,ndm) - h*Cp;
          A(mm+2,ndm) = A(mm+2,ndm) - h*Cq;
          A(mm+3,ndm) = A(mm+3,ndm) - h*Cm;
        end   %end do
        %Add arc length constraint terms to  Hessian and residual
        for k = 1 : ndm 
          b(ndm) = b(ndm) + (x(k) - xo(k))^2;  
          A(ndm,k) = 2*(x(k) - xo(k));
        end   
        b(ndm) = b(ndm) - alpha^2; 
        %norm of residual for convergence test
        test = 0;
        for k = 1 : ndm  
           test = test + b(k)^2;
        end   
        test = sqrt(test);    
        %Update state vector
        x=x-A\b;
        %Test for convergence, if successful output values
        if ((test > tol) & (nu < maxit)) 
            goto_one=1; 
        else 
            goto_one=0; 
        end          
    end %while goto_one==1
    if (nu >= maxit) 
        error('ERROR: Iteration limit exceeded.'); 
    end   
    results(x,column_length,nstep,nu,test,ndm,nbasis, npts,nnstep);   %call results(x,column_length,nstep,nu,test,ndm,nbasis)
    %Update values for previous converged state and guess at next state
    temp = xo;
    xo = x;
    x = 2*x - temp;
end %2 continue
error('Error: Maximum number of steps exhausted');  %stop 'Maximum number of steps exhausted'

function [weight_factor] = simpson(m,npoints,dz,column_length)   %Evaluate the weight weight_factor for the current integration point m |
    c = column_length*dz/3;  
    n = mod(m,2);   
    if((m == 1) | (m == npoints)) 
        weight_factor = c;
    elseif (n==0)
        weight_factor = 4*c;   
    else 
        weight_factor = 2*c; 
    end

function [A,b] = fcn(A,b,D,x,z,column_length,weight_factor,ndm,nbasis, dmo,dqo,dpo) %Compute contribution to K and g at current integration point.
    %displacements, derivatives, and current load weight_factor clf 
    du = 0;
    dw = 0;
    dtheta = 0;
    theta  = 0;
    for i = 1 : nbasis 
        [h,dh] = basis(i,z,column_length,nbasis);    
        du = du + x(3*i-2)*dh;
        dw = dw + x(3*i-1)*dh;
        dtheta = dtheta + x(3*i)*dh;
        theta  = theta  + x(3*i)*h;
    end 
    ct = cos(theta); 
    st = sin(theta); 
        
    %axial strain, shear strain, and curvature
    epsilon = dw*st + (1 + du)*ct - 1;
    beta = dw*ct - (1 + du)*st;
    curv = dtheta;

    %axial force, shear force, bending moment, and other forces    
    bend  = D(1)*curv;
    shear = D(2)*beta;
    axial = D(3)*epsilon;
    Hor = axial*ct - shear*st;
    Ver = axial*st + shear*ct;
    Xi  = (1 + du)*Hor + dw*Ver;
    Yi  = dw*Hor - (1 + du)*Ver;
    %components of E(t)DE + G store in matrix G
    G(1,1) = D(2)*st*st + D(3)*ct*ct;
    G(1,2) = ct*st*(D(3)-D(2));
    G(1,3) = D(2)*st*(1+epsilon) + D(3)*ct*beta - Ver;
    G(1,4) = 0;
    G(2,2) = D(2)*ct*ct + D(3)*st*st;
    G(2,3) = D(3)*st*beta - D(2)*ct*(1+epsilon) + Hor;
    G(2,4) = 0;
    G(3,3) = D(2)*(1+epsilon)^2 + D(3)*beta^2 - Xi;   
    G(3,4) = 0;
    G(4,4) = D(1);
    %rest of G by symmetry
    for i = 1 : 3  
        for j = i+1 : 4  
            G(j,i) = G(i,j); 
        end 
    end 
    %stiffness matrix K, store it in matrix A
    for i = 1 : nbasis 
        [hi,dhi] = basis(i,z,column_length,nbasis);   
        for j = 1 : nbasis 
            [hj,dhj]=basis(j,z,column_length,nbasis); 
            %Bi'GBj noting the sparse structure of B
            for k = 1 : 4  
                GB(k,1) = dhj*G(k,1);   
                GB(k,2) = dhj*G(k,2);   
                GB(k,3) =  hj*G(k,3) + dhj*G(k,4);  
            end 
            for k=1:3 
                BGB(1,k) = dhi*GB(1,k);
                BGB(2,k) = dhi*GB(2,k);
                BGB(3,k) =  hi*GB(3,k) + dhi*GB(4,k);
            end 
            %Assemble the result into the matrix A 
            for m=1:3 
                mm = 3*(i-1);
                for n=1:3    
                  nn = 3*(j-1);
                  A(mm+m,nn+n) = A(mm+m,nn+n)+ BGB(m,n)*weight_factor;
                end   
            end     
        end   
        %Form integral part of residual force and assemble into matrix b
        mm = 3*(i-1);
        [dm,dq,dp]=applied(z,x(ndm),1,dmo,dqo,dpo);  
        b(mm+1) = b(mm+1) + (dhi*Hor  - hi*dp)*weight_factor;
        b(mm+2) = b(mm+2) + (dhi*Ver  - hi*dq)*weight_factor;
        b(mm+3) = b(mm+3) + (dhi*bend + hi*(Yi - dm))*weight_factor;
        %Form the integral part of load weight_factor part of  matrix A
        [dm,dq,dp]=applied(z,x(ndm),2,dmo,dqo,dpo);    
        A(mm+1,ndm) = A(mm+1,ndm) - hi*dp*weight_factor;
        A(mm+2,ndm) = A(mm+2,ndm) - hi*dq*weight_factor;
        A(mm+3,ndm) = A(mm+3,ndm) - hi*dm*weight_factor;
    end %for i 
      
function [h,dh] = basis(i,z,column_length,nbasis) %Evaluate ith basis function h and derivative dh at point z  
    n = mod(i,2);
    m = i/2 + n;
    a = (2*m-1)*2*atan(1); 
    if (n == 0) 
        h = sin(a*z); 
        dh = a*cos(a*z)/column_length;    
    else
        h = 1 - cos(a*z);   
        dh = a*sin(a*z)/column_length;    
    end 

function [dm,dq,dp] = applied(z,clf,n,dmo,dqo,dpo)  %Evaluate the distributed load functions at point z	 
    %nominal values of the applied forces at point z
    f = 1 - z; 
    switch n            
        case 1  %Compute total transverse loads
            dm = dmo(1) + (dmo(2) + dmo(3)*f)*clf;    
            dq = dqo(1) + (dqo(2) + dqo(3)*f)*clf;
            dp = dpo(1) + (dpo(2) + dpo(3)*f)*clf;        
        case 2  %transverse loads associated with load weight_factor only
            dm = dmo(2) + dmo(3)*f;   
            dq = dqo(2) + dqo(3)*f;
            dp = dpo(2) + dpo(3)*f;
    end

function results(x,column_length,nstep,nu,test,ndm,nbasis,  npts,nnstep)  %subroutine results(x,column_length,nstep,nu,test,ndm,nbasis)
    %Determine if current step is an output step
    if(mod(nstep,nnstep) == 0)    
        plot = 1;  
        nplot = npts;
    else
        plot = 0;  
        nplot = 1;
    end 
    %Write coefficients for current step
    if (plot)      
        fprintf('load step=nstep=%g \t load weight_factor=x(ndm)=%g \t iteration=nu=%g \t norm of residual=test=%g \n',nstep, x(ndm), nu, test);  
        for i=1:nbasis   
            ii = 3*(i-1);
            fprintf(' i=%f \t x(ii+1)=%g \t x(ii+2)=%g \t x(ii+3)=%',i,x(ii+1),x(ii+2),x(ii+3));  
        end 
    end 
    %Compute and print current geometry of beam
    if (plot) 
        fprintf('nstep=%g \t x(ndm)=%g\n',nstep,x(ndm)); 
    end
    z = 0;
    dz = 1/nplot;   
    for ii=1:nplot+1
        for j=1:3 
            disp(j) = 0;    
            for i=1:nbasis 
                [h,dh]=basis(i,z,column_length,nbasis); 
                disp(j) = disp(j) + h*x(3*(i-1)+j); 
            end   
        end 
        if (plot) 
            fprintf('z*column_length=%g \t z*column_length+disp(1)=%g \t disp(2)=%g\n',z*column_length,z*column_length+disp(1),disp(2)); 
        end   
        z = z + dz;
    end 
    fprintf('nstep=%g \t x(ndm)=%g \t disp(1)=%g \t disp(2)=%g \t disp(3)=%g \t nu=%g \t test=%g\n',nstep,x(ndm),disp(1),disp(2),disp(3),nu,test);    


    
    