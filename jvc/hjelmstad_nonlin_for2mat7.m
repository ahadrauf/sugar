%hjelmstad_nonlin_for2mat7 

%input:
%    convergence_tolerance=0.1e-7; 
%    alpha=10;
%    column_length=10;
%    num_load_steps=10;
%    max_num_iterations=40;
%    num_integration_intervals=10;
%    num_output_points=20;
%    output_interval=10;
%    num_basis_functions=6;  
%    ndm=3*num_basis_functions+1; 
%    EI=1000; %bend
%    GA=0.1e7; %shear
%    EA=0.1e7; %axial
%    D=[EI;GA;EA]; 
%    Cm=6.28;
%    Cq=0;
%    Cp=0;
%    dmo=[0;0;0];
%    dqo=[0;0;0];
%    dpo=[0;0;0];
%output:
%test1
%Step 	Load 			u(L) 			w(L) 			theta(L) 		NU 	||b||
%1 		9.39690e+000 	-5.66307e-001 	2.85468e+000 	5.85283e-001 	8 	3.84624e-009
%2 		1.88701e+001 	-2.12812e+000 	5.21364e+000 	1.16045e+000 	6 	7.16635e-010
%3 		2.83547e+001 	-4.32027e+000 	6.71375e+000 	1.73261e+000 	7 	1.16761e-009
%4 		3.75724e+001 	-6.89252e+000 	7.20491e+000 	2.35696e+000 	8 	1.33523e-009
%5 		4.70653e+001 	-9.31572e+000 	6.66909e+000 	2.97109e+000 	10 	3.52472e-009
%6 		5.65482e+001 	-1.11767e+001 	5.33214e+000 	3.59404e+000 	8 	2.11990e-009
%7 		6.61260e+001 	-1.20473e+001 	3.66712e+000 	4.15766e+000 	9 	2.27465e-009
%8 		7.58992e+001 	-1.20120e+001 	2.00953e+000 	4.69870e+000 	10 	4.12691e-009
%9 		8.57137e+001 	-1.11950e+001 	5.97090e-001 	5.33377e+000 	8 	2.41631e-009
%10		9.55443e+001 	-1.00686e+001 	2.95210e-003 	5.99692e+000 	10 	2.20560e-009
%time =19.5780

%test 2
%Step 	    Load 			u(L) 			w(L) 			theta(L) 		NU 	||b||
%10 		4.91611e+000 	-1.36094e-002 	4.78034e-001 	6.88999e-002 	5 	9.04853e-010
%20 		9.27598e+000 	-3.72059e-001 	2.43613e+000 	3.82621e-001 	7 	4.29332e-009
%30 		1.15353e+001 	-2.92599e+000 	6.23469e+000 	1.11054e+000 	6 	2.13631e-009
%40 		1.51241e+001 	-6.18260e+000 	7.85273e+000 	1.66807e+000 	6 	1.38255e-009
%50 		1.96751e+001 	-8.48403e+000 	7.98598e+000 	2.00498e+000 	5 	8.93619e-010
%60 		2.44502e+001 	-1.00591e+001 	7.60374e+000 	2.24205e+000 	5 	3.88756e-009
%70 		2.93129e+001 	-1.11290e+001 	7.10690e+000 	2.41766e+000 	5 	3.84823e-009
%80 		3.42428e+001 	-1.18565e+001 	6.64683e+000 	2.54468e+000 	5 	4.85048e-009
%90 		3.92070e+001 	-1.23707e+001 	6.25326e+000 	2.63684e+000 	4 	5.83943e-009
%100 		4.41869e+001 	-1.27503e+001 	5.92188e+000 	2.70459e+000 	4 	2.04446e-009
%120 		5.41644e+001 	-1.32719e+001 	5.40624e+000 	2.79184e+000 	4 	1.01262e-009
%time =192.8570 

function NonlinearBeam 
tic
%   Compute equilibrium path  for Fully Nonlinear Beam      
%      x(3i-2) = a(i), x(3i-1) = b(i), x(3i) = c(i)         
%      a(i) = u coeff, b(i) = w coeff, c(i) = theta coeff   
%      x(ndm) = load factor
%      xo(i) = values of x(i) at last converged state       
%      A = Coefficient matrix, b = right hand side          
%      D = Constitutive properties                          
%      Cm,Cq,Cp Concentrated end forces (M,Q,P)             
%      dmo,dqo,dpo Coefficients for distributed forces      

%Read problem parameters
test0=0;
test1=0;
test2=0;
test3=1;
if test0
    convergence_tolerance=input('convergence convergence_tolerance = ');    
    alpha=input('arc length parameter = ');
    column_length=input('column length = ');
    num_load_steps=input('number of load steps = ');  
    max_num_iterations=input('maximum number of iterations = ');
    num_integration_intervals=input('number of integration intervals = ');
    num_output_points=input('number of points in plotted shape = ');  
    output_interval=input('output shape every N points, N = ');
    D(1)=input('EI = ');  
    D(2)=input('GA = ');
    D(3)=input('EA = ');
    Cm=input('ML = M0 = ');  
    Cq=input('VL = Q0 = ');
    Cp=input('HL = P0 = ');
    dmo=input('distributed load [m01 m02 m03] = dmo(1,2,3) = '); 
    dqo=input('distributed load [q01 q02 q03] = dqo(1,2,3) = '); 
    dpo=input('distributed load [p01 p02 p03] = dpo(1,2,3) = '); 
elseif test1 
    'Pure moment'
    convergence_tolerance=0.1e-7; 
    alpha=10;
    column_length=10;
    num_load_steps=10;
    max_num_iterations=40;
    num_integration_intervals=10;
    num_output_points=20;
    output_interval=10;
    num_basis_functions=6;  
    ndm=3*num_basis_functions+1; 
    EI=1000; %bend
    GA=0.1e7; %shear
    EA=0.1e7; %axial
    D=[EI;GA;EA]; 
    Cm=6.28;
    Cq=0;
    Cp=0;
    dmo=[0;0;0];
    dqo=[0;0;0];
    dpo=[0;0;0];
elseif test2 
    'transverse force and compresive end loads'
    convergence_tolerance=0.1e-7; 
    alpha=0.5;%10;
    column_length=10;
    num_load_steps=120;%10;
    max_num_iterations=40;
    num_integration_intervals=10;
    num_output_points=20;
    output_interval=10;
    num_basis_functions=6;  
    ndm=3*num_basis_functions+1; 
    EI=1000; %bend
    GA=0.1e7; %shear
    EA=0.1e7; %axial
    D=[EI;GA;EA]; 
    Cm=0;%6.28;
    Cq=0;
    Cp=-2.46;%0;
    dmo=[0;0;0];
    dqo=[0.2;0;0];
    dpo=[0;0;0];
elseif test3 
    'Moment and lateral force'
    convergence_tolerance=0.1e-7; 
    alpha=10;
    column_length=10;
    num_load_steps=10;
    max_num_iterations=40;
    num_integration_intervals=10;
    num_output_points=20;
    output_interval=10;
    num_basis_functions=12;  
    ndm=3*num_basis_functions+1; 
    EI=1000; %bend
    GA=0.1e7; %shear
    EA=0.1e7; %axial
    D=[EI;GA;EA]; 
    Cm=0.05;
    Cq=0.05;
    Cp=0.0;
    dmo=[0;0;0];
    dqo=[0;0;0];
    dpo=[0;0;0];
end

fprintf('Step \tLoad \t\t\tu(L) \t\t\tw(L) \t\t\ttheta(L) \t\tNU \t||b||\n'); 
                               
%Initialize 
    dz = 1/(2*num_integration_intervals);         %integration increment
    npoints = 2*num_integration_intervals + 1;    %change intervals to Simpson points
    [x]=zeros(ndm,1);         %load step
    [xo]=zeros(ndm,1);        %load step
    x(ndm) = alpha;           %first load level
          
%Compute MAXSTEPS points along the Equilibrium Path

for nstep=1:num_load_steps   
    nu = 0; %Newton iteration at each load step
    goto_one=1;      
    Cpqm=[Cp;Cq;Cm];
    while goto_one==1      
        nu = nu + 1;        %increment interation counter
        %Numerical integration of Hessian and residual components
        [b]=zeros(ndm,1);   %residual
        [A]=zeros(ndm,ndm); %Hessian
        z = 0;              %Simpson coord
        for mpoint=1:npoints   
          [weight_factor] = simpson(mpoint,npoints,dz,column_length); 
          [A,b] = fcn(A,b,D,x,z,column_length,weight_factor,ndm,num_basis_functions,dmo,dqo,dpo);    
          z = z + dz; %z = z + dz
        end %end do
        %Add end load terms to the residual and coefficient matrix
        for i=1:num_basis_functions   
          [h,dh]=basis(i,1,column_length,num_basis_functions);    
          mm = 3*(i-1);
          b(mm+1:mm+3)= b(mm+1:mm+3) - h*x(ndm)*Cpqm;
          A(mm+1:mm+3,ndm) = A(mm+1:mm+3,ndm) - h*Cpqm;
        end   %end do
        b(ndm) = b(ndm) + sum((x - xo).^2);  
        A(ndm,:) = 2*(x - xo)';       
        b(ndm) = b(ndm) - alpha^2; 
        %norm of residual for convergence test
        test = 0;
        test = test + sum(b.^2);
        test = sqrt(test);    
        %Update state vector
        x=x-A\b;
        %Test for convergence, if successful output values
        if ((test > convergence_tolerance) & (nu < max_num_iterations)) 
            goto_one=1; 
        else 
            goto_one=0; 
        end          
    end %while goto_one==1
%    if (nu >= max_num_iterations) 
%        error('ERROR: Iteration limit exceeded.'); 
%    end   
    results(x,column_length,nstep,nu,test,ndm,num_basis_functions, num_output_points,output_interval);   
    %Update values for previous converged state and guess at next state
    temp = xo;
    xo = x;
    x = 2*x - temp;
end %2 continue
time=toc

function [weight_factor] = simpson(m,npoints,dz,column_length)   %Evaluate the weight weight_factor for the current integration point m 
    %See Hjestad page 386.
    %There are M+1 points in the integration interval, which go from m=0 to m=M.
    M = npoints;    
    c = column_length*dz/3;  
    n = mod(m,2);  %0=even, 1=odd
    if((m == 1) | (m == M)) 
        weight_factor = c; 
    elseif (n==0) %m=odd
        weight_factor = 4*c;   
    else %m=even
        weight_factor = 2*c; 
    end

function [A,b] = fcn(A,b,D,x,z,column_length,weight_factor,ndm,num_basis_functions, dmo,dqo,dpo) %Compute contribution to K and g at current integration point.
    %displacements, derivatives, and current load weight_factor clf 
    du = 0;
    dw = 0;
    dtheta = 0;
    theta  = 0;
    i = 1 : num_basis_functions;
    [h,dh] = Xbasis(i,z,column_length,num_basis_functions);    
    du = du + sum(x(3*i-2)'*dh);
    dw = dw + sum(x(3*i-1)'*dh);
    dtheta = dtheta + sum(x(3*i)'*dh);
    theta  = theta  + sum(x(3*i)'*h);
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
    G=G+G'-diag(diag(G));
    %stiffness matrix K, store it in matrix A
    for i = 1 : num_basis_functions 
        [hi,dhi] = basis(i,z,column_length,num_basis_functions);   
        for j = 1 : num_basis_functions 
            [hj,dhj]=basis(j,z,column_length,num_basis_functions); 
            %Bi'GBj noting the sparse structure of B
            GB(:,1) = dhj*G(:,1);   
            GB(:,2) = dhj*G(:,2);   
            GB(:,3) =  hj*G(:,3) + dhj*G(:,4);  
            BGB(1,:) = dhi*GB(1,:);
            BGB(2,:) = dhi*GB(2,:);
            BGB(3,:) =  hi*GB(3,:) + dhi*GB(4,:);
            %Assemble the result into the matrix A 
            nn = 3*(j-1);
            mm = 3*(i-1);
            A(mm+1:mm+3,nn+1:nn+3) = A(mm+1:mm+3,nn+1:nn+3)+ BGB*weight_factor;
        end   
        %Form integral part of residual force and assemble into matrix b
        mm = 3*(i-1);
        [dm,dq,dp]=applied(z,x(ndm),1,dmo,dqo,dpo);  
        b(mm+1:mm+3) = b(mm+1:mm+3) + [(dhi*Hor  - hi*dp);(dhi*Ver  - hi*dq);(dhi*bend + hi*(Yi - dm))]*weight_factor;
        %Form the integral part of load weight_factor part of  matrix A
        [dm,dq,dp]=applied(z,x(ndm),2,dmo,dqo,dpo);    
        A(mm+1:mm+3,ndm) = A(mm+1:mm+3,ndm) - hi*[dp;dq;dm]*weight_factor;
    end %for i 
      
function [h,dh] = Xbasis(i,z,column_length,num_basis_functions) %Evaluate ith basis function h and derivative dh at point z  
    x1= 1;
    x2= 1.5;
    x3= 1.5;
    for j=1:num_basis_functions
        n = mod(j,2);
        m = floor(j/2) + n;
        a = (2*m-1)*pi/2;
        if (n==0)
            h(j,1) = sin(a*z);
            dh(j,1) = a*cos(a*z)/column_length;
        else
            h(j,1) = 1 - cos(a*z);
            dh(j,1) = a*sin(a*z)/column_length;
        end 
    end
function [h,dh] = basis(i,z,column_length,num_basis_functions) %Evaluate ith basis function h and derivative dh at point z  
    x1= 1;
    x2= 1.5;
    x3= 1.5;
        n = mod(i,2);
        m = floor(i/2) + n;
        a = (2*m-1)*pi/2;
        if (n==0)
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

function results(x,column_length,nstep,nu,test,ndm,num_basis_functions,  num_output_points,output_interval)  
    %Determine if current step is an output step
    if(mod(nstep,output_interval) == 0)    
        plot = 1;  
        nplot = num_output_points;
    else
        plot = 0;  
        nplot = 1;
    end 
    z = 0;
    dz = 1/nplot;   
    for ii=1:nplot+1
        for j=1:3 
            disp(j) = 0;    
%            for i=1:num_basis_functions 
%                [h,dh]=basis(i,z,column_length,num_basis_functions); 
%                disp(j) = disp(j) + h*x(3*(i-1)+j); 
%            end   

            i=1:num_basis_functions;
            [h,dh]=Xbasis(i,z,column_length,num_basis_functions); 
            disp(j) = disp(j) + sum(h'*x(3*(i-1)+j)); 
            
        end 
        z = z + dz;
    end 
    fprintf('%d \t\t%.5e \t%.5e \t%.5e \t%.5e \t%d \t%.5e\n',nstep,x(ndm),disp(1),disp(2),disp(3),nu,test);    

