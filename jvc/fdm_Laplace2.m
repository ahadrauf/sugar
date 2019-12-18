%fdm_Laplace2
%This program solves Laplace's equation u_xx + u_yy = 0 = f(x,y) in the domain xmin < x < xmax,  ymin < y < ymax.
%It has Dirichlet boundary conditions along y=ymin, y=ymax, and a Neumann boundary conditions along x=min and x=max. 
%The sources and boundary conditions are set in the functions at the end of this code. The spacial limits are set by xmin,xmax,ymin,ymax.

function [X,Y,u]=fdm_Laplace2
%By jvclark - summer 2003. 
%From studies of numerical solution of pde's - finite difference methods by Z. Li.

%set boundary limits
gap = 1e-6; 
length = 20e-6;
V0 = 15;
b = 1/50;
a = 100;
param=[gap,V0,b,a,length];
xmin=-length/2; xmax=length/2;    %xmin < x < xmax
ymin=0; ymax=gap;    %ymin < y < ymax

%set the mesh resolution
xres=40*10;   %number of elements along x
yres=80;   %number of elements along y

%set uniform grid size. 
hx = (xmax-xmin)/xres;   hxhx = hx*hx;   %x mesh
hy = (ymax-ymin)/yres;   hyhy = hy*hy;   %y mesh

%create the x and y gridpoints
x=zeros(xres+1,1); y=zeros(yres+1,1);
for i = 1 : xres+1
  x(i) = xmin + (i-1)*hx; %gridpoints along the x-axis
end
for j = 1 : yres+1
  y(j) = ymin + (j-1)*hy; %gridpoints along the y-axis
end

%create the coefficient matrix A, and source matrix f 
M = (yres-1)*xres; %matrix size 
A = sparse(M,M); f = zeros(M,1);
for j = 1 : yres-1 %for each row
    for i = 1 : xres %for each col
        k = i + (j-1)*xres; %matrix-to-vector index
        f(k) = f_source(x(i),y(j+1),param); %evaluate the functions at (x,y)
        A(k,k) = -2/hxhx + -2/hyhy; %coefficient 
        if i == 1 %if x=xmin, at the left x boundary, Neumann BC 
            A(k,k+1) = 2/hxhx;  %coefficient 
            f(k) = f(k) + 2*u_Neumann(y(j+1),param)/hx;    %source
        elseif i==xres %if x=xmax, at the right x boundary, Neumann BC
            A(k,k-1) = 2/hxhx;    %coefficient
            f(k) = f(k) + 2*u_Neumann(y(j+1),param)/hx;    %source
        else
            A(k,k-1) = 1/hxhx;  %coefficient 
            A(k,k+1) = 1/hxhx;  %coefficient 
        end %if x boundary
        if j == 1 %if y=ymin, at the lower y boundary, Dirichlet BC
            A(k,k+xres) = 1/hyhy;
            f(k) = f(k) - u_Dirichlet(x(i),ymin,param)/hyhy;    %source
        elseif j == yres-1 %if y=ymax, at the upper y boundary, Dirichlet BC 
            A(k,k-xres) = 1/hyhy;   %coefficient 
            f(k) = f(k) - u_Dirichlet(x(i),ymax,param)/hyhy;    %source
        else %then not at the boundary, somewhere inside at (x,y) 
            A(k,k-xres) = 1/hyhy;   %coefficient 
            A(k,k+xres) = 1/hyhy;   %coefficient 
        end %if y boundary
    end %for i
end %for j

%solve for the potential inside the region
U = A \ f;

%transform the vector-(i*j,1) form to the matrix-(i,j) form 
j = 1;
for k=1:M
    i = k - (j-1)*xres;
    u(i,j) = U(k);
    u2(i,j) = u_Dirichlet(x(i),y(j+1),param);
    j = fix(k/xres) + 1;
end

%Display results
max_err = max(max(abs(u-u2)))   %the maximum error
X=x(1:xres);  
Y=y(2:yres); 
%plot the solution
figure(1); mesh(Y,X,u); 
title('Solution, U(x,y)'); xlabel('y [meters]'); ylabel('x [meters]'); zlabel('potential [voltage]'); 
%plot the contour
figure(2);contour3(Y,X,u,20)
title('Level curves of U(x,y)'); xlabel('y [meters]'); ylabel('x [meters]'); zlabel('potential [voltage]'); 
figure(3);contour(Y,X,u,20)
title('Level curves of U(x,y)'); xlabel('y [meters]'); ylabel('x [meters]'); zlabel('potential [voltage]'); 

%sources and boundary conditions 
function f = f_source(x,y,param)  %source
    %f = exp(-x)*sin(pi*y)*(1-pi*pi);
    %f = 4;
    f = 0;
    
function f = u_Neumann(y,param)   %boundary condition
    %f = 2;
    %f = -exp(-1)*sin(pi*y);
    f = 0; %dU/dx at y is 0
    
function f = u_Dirichlet(x,y,param) %boundary condition
    %f = x*x + y*y;
    %f = exp(-x)*sin(pi*y);
    gap = param(1);
    V0 = param(2);
    b = param(3);
    a = param(4);
    L = param(5);
    f = y/gap*V0 * (1 + b*sin(a*pi*x/L));


