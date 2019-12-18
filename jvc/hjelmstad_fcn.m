%hjelmstad_fcn 
%function [A,b] = fcn(A,b,D,x,z,column_length,weight_factor,ndm,num_basis_functions, dmo,dqo,dpo) %Compute contribution to K and g at current integration point.

dmo=[0;0;0];
dqo=[0;0;0];
dpo=[0;0;0];
A1=[51324.1551 0 12.3110237 0]';
A2=[0 51324.1551 -136325.467 0]';
A3=[12.3110237 -136325.467 393308.691 0]';
A4=[0 0 0 0]';
A=0*[A1 A2 A3 A4];
b=0*[0.0016019245 18.6549173 -41.4295268 0]';
D(1)=1000;
D(2)=1e6; 
D(3)=1e6;
x=[0 0 0 10]';
z=0.05;
column_length=10;
weight_factor=0.6666666666666666666;
ndm=4;
num_basis_functions=1;


      D
      x
      z
      column_length
      weight_factor
      ndm
      num_basis_functions



%displacements, derivatives, and current load weight_factor clf 
    du = 0;
    dw = 0;
    dtheta = 0;
    theta  = 0;
    for i = 1 : num_basis_functions 

        [h,dh] = hjelmstad_basis(i,z,column_length,num_basis_functions);    
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
    for i = 1 : num_basis_functions 
        [hi,dhi] = hjelmstad_basis(i,z,column_length,num_basis_functions);   
        for j = 1 : num_basis_functions 
            [hj,dhj] = hjelmstad_basis(j,z,column_length,num_basis_functions); 
            %Bi'GBj noting the sparse structure of B
            for k = 1 : 4  
                GB(k,1) = dhj*G(k,1);   
                GB(k,2) = dhj*G(k,2);   
                GB(k,3) =  hj*G(k,3) + dhj*G(k,4);  
            end 
GB(1,1)            
            for k=1:3 
                BGB(1,k) = dhi*GB(1,k);
                BGB(2,k) = dhi*GB(2,k);
                BGB(3,k) =  hi*GB(3,k) + dhi*GB(4,k);
            end 
BGB(1,1)            
            %Assemble the result into the matrix A 
            for m=1:3 
                mm = 3*(i-1);
                for n=1:3    
                  nn = 3*(j-1);
                  mm
                  nn
                  m
                  n
                  weight_factor
                  
a=A(1,1)
bgb=BGB(1,1)
abgb=A(1,1)+BGB(1,1)*weight_factor

                  A(mm+m,nn+n) = A(mm+m,nn+n)+ BGB(m,n)*weight_factor;

A(1,1)
pause        
                  
                end   
            end     
        end   
        
A(1,1)        
pause        
        %Form integral part of residual force and assemble into matrix b
        mm = 3*(i-1);
        [dm,dq,dp] = hjelmstad_applied(z,x(ndm),1,dmo,dqo,dpo);  
        b(mm+1) = b(mm+1) + (dhi*Hor  - hi*dp)*weight_factor;
        b(mm+2) = b(mm+2) + (dhi*Ver  - hi*dq)*weight_factor;
        b(mm+3) = b(mm+3) + (dhi*bend + hi*(Yi - dm))*weight_factor;
        %Form the integral part of load weight_factor part of  matrix A
        [dm,dq,dp] = hjelmstad_applied(z,x(ndm),2,dmo,dqo,dpo);    
        A(mm+1,ndm) = A(mm+1,ndm) - hi*dp*weight_factor;
        A(mm+2,ndm) = A(mm+2,ndm) - hi*dq*weight_factor;
        A(mm+3,ndm) = A(mm+3,ndm) - hi*dm*weight_factor;
    end %for i 
    
A
b

%A =[          57853.6602308377                         0          15.9306003838009  0
%                         0          57853.6602308377         -161798.443869704      0
%          15.9306003838009         -161798.443869704          492690.679711246      0
%                         0                         0                         0      0];
%b =[        0.0020729079700393
%          18.9125654866159
%         -42.4329672106262
%                         0];