%hjelmstad_nonlin_for2mat 

function NonlinearBeam %program NonlinearBeam

%C   *-----------------------------------------------------------*
%C   |   Compute equilibrium path  for Fully Nonlinear Beam      |
%C   |      x(3i-2) = a(i), x(3i-1) = b(i), x(3i) = c(i)         |
%C   |      a(i) = u coeff, b(i) = w coeff, c(i) = theta coeff   |
%C   |      x(ndm) = load factor_                                 |
%C   |      xo(i) = values of x(i) at last converged state       |
%C   |      A = Coefficient matrix, b = right hand side          |
%C   |      D = Constitutive properties                          |
%C   |      Cm,Cq,Cp Concentrated end forces (M,Q,P)             |
%C   |      dmo,dqo,dpo Coefficients for distributed forces      |
%C   *-----------------------------------------------------------*
%      implicit double precision (a-h,o-z)

    nbasis=6;       %parameter (nbasis=6) 
    ndm=3*nbasis+1; %parameter (ndm=3*nbasis+1) 

    A=zeros(ndm,ndm); b=zeros(ndm,1); x=zeros(ndm,1); xo=zeros(ndm,1); D=zeros(3,1); %dimension A(ndm,ndm),b(ndm),x(ndm),xo(ndm),D(3) 

    %common /loads/ dmo(3,1), dqo(3,1), dpo(3,1) 
    %common /out/ npts, nnstep 

    %data zero/0.d0/,one/1.d0/,two/2.d0/

%c.... Read problem parameters
if 0
    tol=input('convergence tolerance = ');    %read(5,1000) tol,alpha,xlength
    alpha=input('arc length parameter = ');
    xlength=input('column length = ');
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
    xlength=10;
    maxsteps=10;
    maxit=40;
    inter=10;
    npts=20;
    nnstep=10;
    D(1)=1000;
    D(2)=0.1e7;
    D(2)=0.1e7;
    Cm=6.28;
    Cq=0;
    Cp=0;
    dmo=[0,0,0];
    dqo=[0,0,0];
    dpo=[0,0,0];
end
%    fprintf('tol=%g \t alpha=%g \t maxsteps=%g \t maxit=%g \t inter=%g \t npts=%g \t nnstep=%g \t nbasis=%g',tol,alpha,maxsteps,maxit,inter,npts,nnstep,nbasis);    %write(6,2000) tol,alpha,maxsteps,maxit,inter,npts,nnstep,nbasis
%    fprintf('xlength=%g \t D(1)=%g \t D(2)=%g \t D(3)=%g \t Cm=%g \t Cq=%g \t Cp=%g',xlength,D(1),D(2),D(3),Cm,Cq,Cp);
%    fprintf('dmo(1)=%g \t dqo(1)=%g \t dpo(1)=%g',dmo(1),dqo(1),dpo(1)); %write(6,2001) xlength,D(1),D(2),D(3),Cm,Cq,Cp,
%    fprintf('dmo(2)=%g \t dqo(2)=%g \t dpo(2)=%g',dmo(2),dqo(2),dpo(2)); % * (dmo(i),dqo(i),dpo(i),i=1,3)
%    fprintf('dmo(3)=%g \t dqo(3)=%g \t dpo(3)=%g',dmo(3),dqo(3),dpo(3))

%c.... Open output files and write header for summary output
%      open(3,file='disp.beam',status='unknown')
%      open(4,file='coef.beam',status='unknown')
    fprintf(' Step \t Load \t u(L) \t w(L) \t theta(L) \t   *  \t NU \t || b ||'); %      write(6,2002)
                               
%c.... Initialize integration increment, change intervals to Simpson points
      dz = 1/(2*inter);   %dz = one/(2*inter)
      npoints = 2*inter + 1;    %npoints = 2*inter + 1

%c.... Initialize values for load step zero,set first load level
      [x]=zeros(ndm,1);   %call zerovec(x,ndm)
      [xo]=zeros(ndm,1);     %call zerovec(xo,ndm)
      x(ndm) = alpha;   %x(ndm) = alpha
          
%c.... Compute MAXSTEPS points along the Equilibrium Path
for nstep=1:maxsteps   %do 2 nstep = 1,maxsteps

%c.... Perform Newton iteration at each load step
      nu = 0; %nu = 0
goto_one=1;      
while goto_one==1      
      nu = nu + 1;    %1 nu = nu + 1
 
%c.... Execute numerical integration of Hessian and residual components
        [b]=zeros(ndm,1);    %call zerovec(b,ndm)
        [A]=zeros(ndm,ndm);    %call zerovec(A,ndm*ndm)
        z = 0;   %z = 0.d0
        for mpoint=1:npoints   %do mpoint = 1,npoints
          [factor_]=simpson(mpoint,npoints,dz,xlength); %call simpson(mpoint,npoints,dz,xlength,factor_)
          [A,b,D,x,z,xlength,factor_,ndm,nbasis]=fcn(A,b,D,x,z,xlength,factor_,ndm,nbasis,dmo,dqo,dpo);    %call fcn(A,b,D,x,z,xlength,factor_,ndm,nbasis)
          z = z + dz; %z = z + dz
        end %end do

%C.... Add end load terms to the residual and coefficient matrix
        for i=1:nbasis   %do i=1,nbasis
          [h,dh]=basis(i,1,xlength,nbasis);    %call basis(i,one,xlength,h,dh,nbasis)
          mm = 3*(i-1);
          b(mm+1) = b(mm+1) - h*Cp*x(ndm);
          b(mm+2) = b(mm+2) - h*Cq*x(ndm);
          b(mm+3) = b(mm+3) - h*Cm*x(ndm);
          A(mm+1,ndm) = A(mm+1,ndm) - h*Cp;
          A(mm+2,ndm) = A(mm+2,ndm) - h*Cq;
          A(mm+3,ndm) = A(mm+3,ndm) - h*Cm;
        end   %end do

%c.... Add arc length constraint terms to  Hessian and residual
        for k=1:ndm %do k=1,ndm
          b(ndm) = b(ndm) + (x(k) - xo(k))^2;  %b(ndm) = b(ndm) + (x(k) - xo(k))**2
          A(ndm,k) = 2*(x(k) - xo(k));
        end   %end do
        b(ndm) = b(ndm) - alpha^2; %b(ndm) = b(ndm) - alpha**2
 
%c.... Compute norm of residual for convergence test
        test = 0;
        for k=1:ndm  %do k=1,ndm
           test = test + b(k)^2;   %test = test + b(k)**2
        end   %end do
        test = sqrt(test);     %test = dsqrt(test)

%c.... Compute eigenvalues of Tangent stiffness matrix
%c        call eigen(A,eval,evec,ndm)

%c.... Update state vector

        [A]=inv(A);  %invert(A,ndm,ndm)
        for j=1:ndm  %do j=1,ndm
          for k=1:ndm    %do k=1,ndm
            x(j) = x(j) - A(j,k)*b(k);
          end %end do
        end %end do

%c.... Test for convergence, if successful output values
        if ((test > tol) & (nu < maxit)) goto_one=1; else goto_one=0; end  %if ((test.gt.tol).and.(nu.lt.maxit)) go to 1
end %while goto_one==1      
        if (nu >= maxit) error('ERROR: Iteration limit exceeded.'); end   %if (nu.ge.maxit) stop 'Iteration limit exceeded'
        results(x,xlength,nstep,nu,test,ndm,nbasis, npts,nnstep);   %call results(x,xlength,nstep,nu,test,ndm,nbasis)
   
%c.... Update values for previous converged state and guess at next state
        for j=1:ndm  %do j=1,ndm
          temp  = xo(j);
          xo(j) = x(j);
          x(j)  = 2*x(j) - temp;
        end %end do

end %2 continue
%      close(3)
%      close(4)
      error('Error: Maximum number of steps exhausted');  %stop 'Maximum number of steps exhausted'

% 1000 format(3f10.0)
% 1001 format(3i10)                                        
% 2000 format(/' Fully Nonlinear Beam Analysis'//
%     *     5x,  ' Convergence tolerance             : ',e12.4/
%     *     5x,  ' Arc length parameter              : ',f12.4/
%     *     5x,  ' Number of load steps              : ',i6/
%     *     5x,  ' Maximum number of iterations      : ',i6/
%     *     5x,  ' Number of integration intervals   : ',i6/
%     *     5x,  ' Number of points in plotted shape : ',i6/
%     *     5x,  ' Output shape every N points, N =  : ',i6/
%     *     5x,  ' Number of basis functions         : ',i6/)
% 2001 format(/' Beam properties, end load and dist. load amplitudes'/
%     *     5x,  ' Column Length  : ',g12.4/
%     *     5x,  ' EA, GA, EI     : ',3g12.4/
%     *     5x,  ' Mo, Qo, Po     : ',3g12.4/
%     *     5x,  ' mo1, qo1, po1  : ',3g12.4/
%     *     5x,  ' mo2, qo2, po2  : ',3g12.4/
%     *     5x,  ' mo3, qo3, po3  : ',3g12.4//)
% 2002 format(' Step',8x,'Load',8x,'u(L)',8x,'w(L)',4x,'theta(L)',
%     *       '   NU',4x,'|| b ||')
%      end

%C----------------------------------------------------------------------FCN

function [A,b,D,x,z,xlength,factor_,ndm,nbasis]=fcn(A,b,D,x,z,xlength,factor_,ndm,nbasis, dmo,dqo,dpo)   %subroutine fcn(A,b,D,x,z,xlength,factor_,ndm,nbasis)

%C  
%*------------------------------------------------------------------*
%C   |  Compute contribution to K and g at current integration point | 
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)

%      dimension A(ndm,ndm),b(ndm),x(ndm),D(3)   
%      dimension G(4,4),GB(4,3),BGB(3,3)

%      data zero/0.d0/,one/1.d0/
     
%c.... Compute displacements, derivatives, and current load factor_ clf 

      du = 0;
      dw = 0;
      dtheta = 0;
      theta  = 0;
      for i=1:nbasis %do i=1,nbasis
         [h,dh]=basis(i,z,xlength,nbasis);    %call basis(i,z,xlength,h,dh,nbasis) 
         du = du + x(3*i-2)*dh;
         dw = dw + x(3*i-1)*dh;
         dtheta = dtheta + x(3*i)*dh;
         theta  = theta  + x(3*i)*h;
      end    %end do
      ct = cos(theta); %ct = dcos(theta)
      st = sin(theta); %st = dsin(theta)
        
%c.... Compute axial strain, shear strain, and curvature
      epsilon = dw*st + (1 + du)*ct - 1;
      beta    = dw*ct - (1 + du)*st;
      curv    = dtheta;

%c.... Compute axial force, shear force, bending moment, and other forces
      bend  = D(1)*curv;
      shear = D(2)*beta;
      axial = D(3)*epsilon;
      Hor = axial*ct - shear*st;
      Ver = axial*st + shear*ct;
      Xi  = (1 + du)*Hor + dw*Ver;
      Yi  = dw*Hor - (1 + du)*Ver;

%c.... Compute components of E(t)DE + G store in matrix G
      G(1,1) = D(2)*st*st + D(3)*ct*ct;
      G(1,2) = ct*st*(D(3)-D(2));
      G(1,3) = D(2)*st*(1+epsilon) + D(3)*ct*beta - Ver;
      G(1,4) = 0;
      G(2,2) = D(2)*ct*ct + D(3)*st*st;
      G(2,3) = D(3)*st*beta - D(2)*ct*(1+epsilon) + Hor;
      G(2,4) = 0;
      G(3,3) = D(2)*(1+epsilon)^2 + D(3)*beta^2 - Xi;   %G(3,3) = D(2)*(one+epsilon)**2 + D(3)*beta**2 - Xi
      G(3,4) = 0;
      G(4,4) = D(1);

%c.... Compute the rest of G by symmetry
      for i=1:3  %do i=1,3
        for j=i+1:4  %do j=i+1,4
          G(j,i) = G(i,j);
        end %end do
      end %end do

%c.... Form stiffness matrix K and store it in matrix A
      for i=1:nbasis %do i=1,nbasis
        [hi,dhi]=basis(i,z,xlength,nbasis);   %call basis(i,z,xlength,hi,dhi,nbasis)
        for j=1:nbasis   %do j=1,nbasis
          [hj,dhj]=basis(j,z,xlength,nbasis); %call basis(j,z,xlength,hj,dhj,nbasis)
 
%c.... Compute Bi(t)GBj noting the sparse structure of B
          for k=1:4  %do k=1,4
            GB(k,1) = dhj*G(k,1);
            GB(k,2) = dhj*G(k,2);
            GB(k,3) =  hj*G(k,3) + dhj*G(k,4);
          end %end do
          for k=1:3  %do k=1,3
            BGB(1,k) = dhi*GB(1,k);
            BGB(2,k) = dhi*GB(2,k);
            BGB(3,k) =  hi*GB(3,k) + dhi*GB(4,k);
          end %end do

%c.... Assemble the result into the matrix A 
          for m=1:3  %do m=1,3
            mm = 3*(i-1);
            for n=1:3    %do n=1,3
              nn = 3*(j-1);
              A(mm+m,nn+n) = A(mm+m,nn+n)+ BGB(m,n)*factor_;
            end   %end do
          end     %end do 

        end   %end do %i nbasis

%c.... Form integral part of residual force and assemble into matrix b
        mm = 3*(i-1);
        [dm,dq,dp]=applied(z,x(ndm),1, dmo,dqo,dpo);  %call applied(z,dm,dq,dp,x(ndm),1)
        b(mm+1) = b(mm+1) + (dhi*Hor  - hi*dp)*factor_;
        b(mm+2) = b(mm+2) + (dhi*Ver  - hi*dq)*factor_;
        b(mm+3) = b(mm+3) + (dhi*bend + hi*(Yi - dm))*factor_;

%c.... Form the integral part of load factor_ part of  matrix A
        [dm,dq,dp]=applied(z,x(ndm),2, dmo,dqo,dpo);    %call applied(z,dm,dq,dp,x(ndm),2)
        A(mm+1,ndm) = A(mm+1,ndm) - hi*dp*factor_;
        A(mm+2,ndm) = A(mm+2,ndm) - hi*dq*factor_;
        A(mm+3,ndm) = A(mm+3,ndm) - hi*dm*factor_;

      end %end do %j nbasis          
%      return
%      end

%C--------------------------------------------------------------------BASIS

function [h,dh]=basis(i,z,xlength,nbasis)    %subroutine basis(i,z,xlength,h,dh,nbasis) 		            
      
%C  
%*------------------------------------------------------------------*
%C   |  Evaluate ith basis function h and derivative dh at point z  |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)

      n = mod(i,2);
      m = i/2 + n;
      a = (2*m-1)*2*atan(1); %a = dble(2*m-1)*2.d0*datan(1.d0)
      if (n == 0) %if (n.eq.0) then
        h = sin(a*z);  %h = dsin(a*z)
        dh = a*cos(a*z)/xlength;    %dh = a*dcos(a*z)/xlength
      else
        h = 1 - cos(a*z);   %h = 1.d0 - dcos(a*z)
        dh = a*sin(a*z)/xlength;    %dh = a*dsin(a*z)/xlength
      end %end if
%      return
%      end

%C--------------------------------------------------------------------BASIS

function [h,dh]=basis2(i,z,xlength,h,dh,nbasis)    %subroutine basis2(i,z,xlength,h,dh,nbasis)		             

%C  
%*------------------------------------------------------------------*
%C   |  Evaluate ith basis function h and derivative dh at point z  |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)

      if (nbasis > 5) error('Not enough basis functions in BASIS'); end   %if (nbasis.gt.5) Stop 'Not enough basis functions in BASIS'
      %go to (1,2,3,4,5),i
    switch i

%c.... Basis function h1(z) = z
    case 1
      h  = z;   %1 h  = z
      dh = 1/xlength; %dh = 1.d0/xlength
%      return

%c.... Basis function h2(z) = 4z^2 - 3z
    case 2
      h  = (4*z - 3)*z;   %2 h  = (4.d0*z - 3.d0)*z
      dh = (8*z - 3)/xlength; %dh = (8.d0*z - 3.d0)/xlength
%      return

%c.... Basis function h3(z) = 15z^3 - 20z^2 +6z
    case 3
      h  = ((15*z - 20)*z + 6)*z;  %3 h  = ((15.d0*z - 20.d0)*z + 6.d0)*z
      dh = ((45*z - 40)*z + 6)/xlength;    %dh = ((45.d0*z - 40.d0)*z + 6.d0)/xlength
%      return 

%c.... Basis function h4(z) = 56z^4 - 105z^3 + 60z^2 - 10z
    case 4
      h  = ((( 56*z - 105)*z +  60)*z - 10)*z;  %4 h  = ((( 56.d0*z - 105.d0)*z +  60.d0)*z - 10.d0)*z
      dh = (((224*z - 315)*z + 120)*z - 10)/xlength;    %dh = (((224.d0*z - 315.d0)*z + 120.d0)*z - 10.d0)/xlength
%      return

%c.... Basis function h5(z) = z^5
    case 5
      h  = z^5;    %5 h  = z**5
      dh = (5*z^4)/xlength;     %dh = (5.d0*z**4)/xlength
%      return
     
    end
     
%C-------------------------------------------------------------------APPLIED

function [dm,dq,dp]=applied(z,clf,n,  dmo,dqo,dpo)  %subroutine applied(z,dm,dq,dp,clf,n)                            

%C  
%*------------------------------------------------------------------*
%C   |  Evaluate the distributed load functions at point z	  |
%C   |  DM = moment, DQ = transverse force, DP = axial force  |
%C   |  CLF is the current load factor_.				   |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)
%      common /loads/ dmo(3),dqo(3),dpo(3)

%c.... Compute the nominal values of the applied forces at point z
      f = 1 - z; %f = 1.d0 - z
%      go to (1,2), n
    switch n
%c.... Compute total transverse loads
    case 1
      dm = dmo(1) + (dmo(2) + dmo(3)*f)*clf;    %1 dm = dmo(1) + (dmo(2) + dmo(3)*f)*clf
      dq = dqo(1) + (dqo(2) + dqo(3)*f)*clf;
      dp = dpo(1) + (dpo(2) + dpo(3)*f)*clf;
%      return

%c.... Compute transverse loads associated with load factor_ only
    case 2
      dm = dmo(2) + dmo(3)*f;   %2 dm = dmo(2) + dmo(3)*f
      dq = dqo(2) + dqo(3)*f;
      dp = dpo(2) + dpo(3)*f;
%      return
    end

%C-------------------------------------------------------------------RESULTS

function results(x,xlength,nstep,nu,test,ndm,nbasis,  npts,nnstep)  %subroutine results(x,xlength,nstep,nu,test,ndm,nbasis)

%C  
%*------------------------------------------------------------------*
%C   |  Print results of current step to various files	 |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)
%      logical plot
%      dimension x(ndm),disp(3)
%      common /out/ npts,nnstep

%c.... Determine if current step is an output step
      if(mod(nstep,nnstep) == 0)    %if(mod(nstep,nnstep).eq.0) then
        plot = 1;  %plot = .true.
        nplot = npts;
      else
        plot = 0;     %plot = .false.
        nplot = 1;
      end %end if

%c.... Write coefficients for current step
      if (plot)      %if (plot) then
        fprintf('load step=nstep=%g \t load factor_=x(ndm)=%g \t iteration=nu=%g \t norm of residual=test=%g ',nstep, x(ndm), nu, test);  %write(4,2001) nstep, x(ndm), nu, test
        for i=1:nbasis   %do i=1,nbasis
          ii = 3*(i-1);
          fprintf(' i=%f \t x(ii+1)=%g \t x(ii+2)=%g \t x(ii+3)=%',i,x(ii+1),x(ii+2),x(ii+3));  %write(4,2002) i,(x(ii+k),k=1,3);
        end %end do
      end %end if

%c.... Compute and print current geometry of beam
      if (plot) fprintf('nstep=%g \t x(ndm)=%g',nstep,x(ndm)); end  %write(3,3001) nstep,x(ndm)
      z = 0; %z = 0.d0
      dz = 1/nplot;   %dz = 1.0/nplot
      for ii=1:nplot+1   %do ii=1,nplot+1
        for j=1:3    %do j=1,3
          disp(j) = 0;    %disp(j) = 0.d0
          for i=1:nbasis %do i=1,nbasis
            [h,dh]=basis(i,z,xlength,nbasis); %call basis(i,z,xlength,h,dh,nbasis)
            disp(j) = disp(j) + h*x(3*(i-1)+j); %disp(j) = disp(j) + h*x(3*(i-1)+j)
          end   %end do
        end %end do
        if (plot) fprintf('z*xlength=%g \t z*xlength+disp(1)=%g \t disp(2)=%g',z*xlength,z*xlength+disp(1),disp(2)); end   %write(3,3000) z*xlength,z*xlength+disp(1),disp(2)
        z = z + dz;
      end %end do
      fprintf('nstep=%g \t x(ndm)=%g \t disp(1)=%g \t disp(2)=%g \t disp(3)=%g \t nu=%g \t test=%g',nstep,x(ndm),disp(1),disp(2),disp(3),nu,test);    %write(6,2003) nstep,x(ndm),(disp(j),j=1,3),nu,test
%      return

% 2001 format(/' Load Step  :',i4,5x,' Load factor_      :',e12.5,
%     *       /' Iterations :',i4,5x,' Norm of Residual :',e12.5/
%     *       /'    i',15x,' a(i)',15x,' b(i)',15x,' c(i)')
% 2002 format(i5,3e20.5)
% 2003 format(i5,4e12.5,i5,e12.4)
% 3000 format(4e15.5)
% 3001 format(' Step : ',i5,5x,' Load : ',e15.5)
%      end
 
%C-------------------------------------------------------------------SIMPSON

function [factor_] = simpson(m,npoints,dz,xlength)   %subroutine simpson(m,npoints,dz,xlength,factor_)

%C  
%*------------------------------------------------------------------*
%C   |  Evaluate the weight factor_ for the current integration point m |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)

      c = xlength*dz/3;  %c = xlength*dz/3.d0
      n = mod(m,2);   
      if((m == 1) | (m == npoints)) %if((m.eq.1).or.(m.eq.npoints)) then
         factor_ = c;
      elseif (n==0) %else if (n.eq.0) then
         factor_ = 4*c;   %factor_ = 4.d0*c
      else 
         factor_ = 2*c;    %factor_ = 2.d0*c
      end %endif
%      return
%      end

%C------------------------------------------------------------------ZEROVEC

%      function [v]=zerovec(v,n)
%C  
%*------------------------------------------------------------------*
%C   |  Initialize the array V to zero for all values up to n	 |
%C  
%*------------------------------------------------------------------*
%      implicit double precision (a-h,o-z)
%      dimension v(n)

%      for i=1:n  %do i=1,n
%        v(i) = 0;    %v(i) = 0.d0
%      end %end do
%      return
%      end
  
%C--------------------------------------------------------------------INVERT

%      function [a,nmax,ndm]=invert(a,nmax,ndm)
%C   *---------------------------------------------------------------*
%C   |     Invert matrix A (nmax,nmax) when array dimension is NDM   |
%C   *---------------------------------------------------------------*
%      implicit double precision (a-h,o-z)
%      dimension a(ndm,ndm)

%      for n = 1:nmax %do 200 n = 1,nmax
%        d = a(n,n);
%        for j = 1:nmax   %do 100 j = 1,nmax
%          a(n,j) = -a(n,j)/d;
%        end   %100   continue
%        for i = 1:nmax   %do 150 i = 1,nmax
%          if(n~=i) %if(n.eq.i) go to 150
%             for j = 1:nmax %do 140 j = 1,nmax
%               if(n~=j) a(i,j) = a(i,j) + a(i,n)*a(n,j); %if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
%             end  %140     continue
%             a(i,n) = a(i,n)/d;
%          end
%        end   %150   continue
%        a(n,n) = 1/d;  %a(n,n) = 1.d0/d
%      end   %200 continue
%      return
%      end

