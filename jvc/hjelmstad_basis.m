function [h,dh] = basis(i,z,column_length,num_basis_functions) %Evaluate ith basis function h and derivative dh at point z  


%    i=1;
%    z=0.9;
%    column_length=10;
%    num_basis_functions=6;



    n = mod(i,2);
    m = floor(i/2) + n;
    a = (2*m-1)*2*atan(1); 
    if (n == 0) 
        h = sin(a*z); 
        dh = a*cos(a*z)/column_length;    
    else
        h = 1 - cos(a*z);   
        dh = a*sin(a*z)/column_length;    
    end 

%    i,z,column_length,h,dh,num_basis_functions,n,m,a
