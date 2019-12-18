%fdm_poisson2
%This program solves Poisson equation u_xx + u_yy = f(x,y), a < x < b,  c < y < d.
%with Dirichlet BCs along x=b, y=c, y=d, and a Neumann BC at x=a. 

function fdm_poisson2
%By jvclark - summer 2003. 
%From by studies of numerical solution of pde's - finite difference methose, Zhilin Li.

%set boundary limits
a=1;    b=2;    %a<x<b
c=-1;   d=1;    %c<y<d

%set the mesh resolution
m=40;   %along m
n=80;   %along n

%set uniform grid size. 
hx = (b-a)/m;   hxhx = hx*hx;   %x mesh
hy = (d-c)/n;   hyhy = hy*hy;   %y mesh

%create the x and y gridpoints
x=zeros(m+1,1); y=zeros(n+1,1);
for i=1:m+1
  x(i) = a + (i-1)*hx; %gridpoints along the x-axis
end
for i=1:n+1
  y(i) = c + (i-1)*hy; %gridpoints along the y-axis
end

%create the coefficient matrix A, and source matrix f 
M = (n-1)*m; %matrix size 
A = sparse(M,M); f = zeros(M,1);
for j = 1 : n-1 %for each row
    for i = 1 : m %for each col
        k = i + (j-1)*m; %matrix-to-vector index
        f(k) = f_source(x(i),y(j+1)); %evaluate the functions at (x,y)
        A(k,k) = -2/hxhx -2/hyhy; %coefficient 
        if i == 1 %if x=a, at the left x boundary 
            A(k,k+1) = 2/hxhx;
            f(k) = f(k) + 2*ux_BC(y(j+1))/hx;
        elseif i==m %if x=b, at the right x boundary
            A(k,k-1) = 1/hxhx;
            f(k) = f(k) - ue_BC(x(i+1),y(j+1))/hxhx;
        else
            A(k,k-1) = 1/hxhx;
            A(k,k+1) = 1/hxhx;
        end %if x boundary
        if j == 1 %if y=c, at the lower y boundary
            A(k,k+m) = 1/hyhy;
            f(k) = f(k) - ue_BC(x(i),c)/hyhy;
        elseif j == n-1 %if y=d, at the upper y boundary
            A(k,k-m) = 1/hyhy;
            f(k) = f(k) - ue_BC(x(i),d)/hyhy;
        else
            A(k,k-m) = 1/hyhy;
            A(k,k+m) = 1/hyhy;
        end %if y boundary
    end %for i
end %for j

%solve for the potential inside the region
U = A \ f;

%transform the vector-(i*j,1) form to the matrix-(i,j) form 
j = 1;
for k=1:M
    i = k - (j-1)*m ;
    u(i,j) = U(k);
    u2(i,j) = ue_BC(x(i),y(j+1));
    j = fix(k/m) + 1;
end

%Display results
max_err = max(max(abs(u-u2)))   %the maximum error
X=x(1:m); 
Y=y(2:n);
%plot the solution
figure(1); mesh(Y,X,u); 
title('The solution plot'); xlabel('y'); ylabel('x');
%plot the error
figure(2); mesh(Y,X,u-u2); 
title('The error plot'); xlabel('y'); ylabel('x');

%sources and boundary conditions 
function f = f_source(x,y)  %source
    %f = exp(-x)*sin(pi*y)*(1-pi*pi);
    %f = 4;
    f = 0;
    
function f = ux_BC(y)   %boundary condition
    %f = 2;
    f = -exp(-1)*sin(pi*y);
    
function f = ue_BC(x,y) %boundary condition
    %f = x*x + y*y;
    f = exp(-x)*sin(pi*y);

    
    
    
    