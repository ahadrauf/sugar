%      program NonlinearBeam
%C   *-----------------------------------------------------------*
%C   |   Compute equilibrium path  for Fully Nonlinear Beam      |
%C   |      x(3i-2) = a(i), x(3i-1) = b(i), x(3i) = c(i)         |
%C   |      a(i) = u coeff, b(i) = w coeff, c(i) = theta coeff   |
%C   |      x(ndm) = load factor                                 |
%C   |      xo(i) = values of x(i) at last converged state       |
%C   |      A = Coefficient matrix, b = right hand side          |
%C   |      D = Constitutive properties                          |
%C   |      Cm,Cq,Cp Concentrated end forces (M,Q,P)             |
%C   |      dmo,dqo,dpo Coefficients for distributed forces      |
%C   *-----------------------------------------------------------*
      implicit double precision (a-h,o-z)

      parameter (nbasis=6)
      parameter (ndm=3*nbasis+1)

      dimension A(ndm,ndm),b(ndm),x(ndm),xo(ndm),D(3)
      common /loads/ dmo(3),dqo(3),dpo(3)
      common /out/   npts,nnstep

      data zero/0.d0/,one/1.d0/,two/2.d0/

c.... Read problem parameters
      read(5,1000) tol,alpha,xlength
      read(5,1001) maxsteps,maxit,inter
      read(5,1001) npts,nnstep
      read(5,1000) D(1),D(2),D(3)
      read(5,1000) Cm,Cq,Cp
      read(5,1000) (dmo(i),i=1,3)
      read(5,1000) (dqo(i),i=1,3)
      read(5,1000) (dpo(i),i=1,3)
      write(6,2000) tol,alpha,maxsteps,maxit,inter,npts,nnstep,nbasis
      write(6,2001) xlength,D(1),D(2),D(3),Cm,Cq,Cp,
     *              (dmo(i),dqo(i),dpo(i),i=1,3)

c.... Open output files and write header for summary output
      open(3,file='disp.beam',status='unknown')
      open(4,file='coef.beam',status='unknown')
      write(6,2002)
                                
c.... Initialize integration increment, change intervals to Simpson
points
      dz = one/(2*inter)  
      npoints = 2*inter + 1

c.... Initialize values for load step zero,set first load level
      call zerovec(x,ndm)
      call zerovec(xo,ndm)
      x(ndm) = alpha
          
c.... Compute MAXSTEPS points along the Equilibrium Path
      do 2 nstep = 1,maxsteps

c.... Perform Newton iteration at each load step
        nu = 0
    1   nu = nu + 1
 
c.... Execute numerical integration of Hessian and residual components
        call zerovec(b,ndm)
        call zerovec(A,ndm*ndm)
        z = 0.d0
        do mpoint = 1,npoints
          call simpson(mpoint,npoints,dz,xlength,factor)
          call fcn(A,b,D,x,z,xlength,factor,ndm,nbasis)
          z = z + dz
        end do

C.... Add end load terms to the residual and coefficient matrix
        do i=1,nbasis
          call basis(i,one,xlength,h,dh,nbasis)
          mm = 3*(i-1)
          b(mm+1) = b(mm+1) - h*Cp*x(ndm)
          b(mm+2) = b(mm+2) - h*Cq*x(ndm)
          b(mm+3) = b(mm+3) - h*Cm*x(ndm)
          A(mm+1,ndm) = A(mm+1,ndm) - h*Cp
          A(mm+2,ndm) = A(mm+2,ndm) - h*Cq
          A(mm+3,ndm) = A(mm+3,ndm) - h*Cm
        end do

c.... Add arc length constraint terms to  Hessian and residual
        do k=1,ndm
          b(ndm)   = b(ndm) + (x(k) - xo(k))**2
          A(ndm,k) = two*(x(k) - xo(k))
        end do
        b(ndm) = b(ndm) - alpha**2
 
c.... Compute norm of residual for convergence test
        test = zero
        do k=1,ndm
           test = test + b(k)**2
        end do
        test = dsqrt(test)

c.... Compute eigenvalues of Tangent stiffness matrix
c        call eigen(A,eval,evec,ndm)

c.... Update state vector
        call invert(A,ndm,ndm)
        do j=1,ndm
          do k=1,ndm
            x(j) = x(j) - A(j,k)*b(k)
          end do
        end do

c.... Test for convergence, if successful output values
        if ((test.gt.tol).and.(nu.lt.maxit)) go to 1
        if (nu.ge.maxit) stop 'Iteration limit exceeded'
        call results(x,xlength,nstep,nu,test,ndm,nbasis)
   
c.... Update values for previous converged state and guess at next
state
        do j=1,ndm
          temp  = xo(j)
          xo(j) = x(j)
          x(j)  = two*x(j) - temp
        end do

    2 continue
      close(3)
      close(4)
      stop 'Maximum number of steps exhausted'

 1000 format(3f10.0)
 1001 format(3i10)                                        
 2000 format(/' Fully Nonlinear Beam Analysis'//
     *     5x,  ' Convergence tolerance             : ',e12.4/
     *     5x,  ' Arc length parameter              : ',f12.4/
     *     5x,  ' Number of load steps              : ',i6/
     *     5x,  ' Maximum number of iterations      : ',i6/
     *     5x,  ' Number of integration intervals   : ',i6/
     *     5x,  ' Number of points in plotted shape : ',i6/
     *     5x,  ' Output shape every N points, N =  : ',i6/
     *     5x,  ' Number of basis functions         : ',i6/)
 2001 format(/' Beam properties, end load and dist. load amplitudes'/
     *     5x,  ' Column Length  : ',g12.4/
     *     5x,  ' EA, GA, EI     : ',3g12.4/
     *     5x,  ' Mo, Qo, Po     : ',3g12.4/
     *     5x,  ' mo1, qo1, po1  : ',3g12.4/
     *     5x,  ' mo2, qo2, po2  : ',3g12.4/
     *     5x,  ' mo3, qo3, po3  : ',3g12.4//)
 2002 format(' Step',8x,'Load',8x,'u(L)',8x,'w(L)',4x,'theta(L)',
     *       '   NU',4x,'|| b ||')
      end

C----------------------------------------------------------------------FCN

      subroutine fcn(A,b,D,x,z,xlength,factor,ndm,nbasis)
C  
*------------------------------------------------------------------*
C   |  Compute contribution to K and g at current integration point  
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)

      dimension A(ndm,ndm),b(ndm),x(ndm),D(3)   
      dimension G(4,4),GB(4,3),BGB(3,3)

      data zero/0.d0/,one/1.d0/
     
c.... Compute displacements, derivatives, and current load factor clf
      du = zero
      dw = zero
      dtheta = zero
      theta  = zero
      do i=1,nbasis
         call basis(i,z,xlength,h,dh,nbasis)
         du = du + x(3*i-2)*dh
         dw = dw + x(3*i-1)*dh
         dtheta = dtheta + x(3*i)*dh
         theta  = theta  + x(3*i)*h
      end do
      ct = dcos(theta)
      st = dsin(theta)
        
c.... Compute axial strain, shear strain, and curvature
      epsilon = dw*st + (one + du)*ct - one
      beta    = dw*ct - (one + du)*st
      curv    = dtheta

c.... Compute axial force, shear force, bending moment, and other
forces
      bend  = D(1)*curv
      shear = D(2)*beta
      axial = D(3)*epsilon
      Hor = axial*ct - shear*st
      Ver = axial*st + shear*ct
      Xi  = (one + du)*Hor + dw*Ver
      Yi  = dw*Hor - (one + du)*Ver 

c.... Compute components of E(t)DE + G store in matrix G
      G(1,1) = D(2)*st*st + D(3)*ct*ct
      G(1,2) = ct*st*(D(3)-D(2))
      G(1,3) = D(2)*st*(one+epsilon) + D(3)*ct*beta - Ver
      G(1,4) = zero
      G(2,2) = D(2)*ct*ct + D(3)*st*st
      G(2,3) = D(3)*st*beta - D(2)*ct*(one+epsilon) + Hor
      G(2,4) = zero
      G(3,3) = D(2)*(one+epsilon)**2 + D(3)*beta**2 - Xi
      G(3,4) = zero
      G(4,4) = D(1)

c.... Compute the rest of G by symmetry
      do i=1,3
        do j=i+1,4
          G(j,i) = G(i,j)
        end do
      end do

c.... Form stiffness matrix K and store it in matrix A
      do i=1,nbasis
        call basis(i,z,xlength,hi,dhi,nbasis)
        do j=1,nbasis
          call basis(j,z,xlength,hj,dhj,nbasis)
 
c.... Compute Bi(t)GBj noting the sparse structure of B
          do k=1,4
            GB(k,1) = dhj*G(k,1)
            GB(k,2) = dhj*G(k,2)
            GB(k,3) =  hj*G(k,3) + dhj*G(k,4)
          end do
          do k=1,3
            BGB(1,k) = dhi*GB(1,k)
            BGB(2,k) = dhi*GB(2,k)
            BGB(3,k) =  hi*GB(3,k) + dhi*GB(4,k)
          end do

c.... Assemble the result into the matrix A 
          do m=1,3
            mm = 3*(i-1)
            do n=1,3
              nn = 3*(j-1)
              A(mm+m,nn+n) = A(mm+m,nn+n)+ BGB(m,n)*factor
            end do
          end do 

        end do

c.... Form integral part of residual force and assemble into matrix b
        mm = 3*(i-1)
        call applied(z,dm,dq,dp,x(ndm),1)
        b(mm+1) = b(mm+1) + (dhi*Hor  - hi*dp)*factor
        b(mm+2) = b(mm+2) + (dhi*Ver  - hi*dq)*factor
        b(mm+3) = b(mm+3) + (dhi*bend + hi*(Yi - dm))*factor

c.... Form the integral part of load factor part of  matrix A
        call applied(z,dm,dq,dp,x(ndm),2)
        A(mm+1,ndm) = A(mm+1,ndm) - hi*dp*factor
        A(mm+2,ndm) = A(mm+2,ndm) - hi*dq*factor
        A(mm+3,ndm) = A(mm+3,ndm) - hi*dm*factor

      end do          
      return
      end

C--------------------------------------------------------------------BASIS

      subroutine basis(i,z,xlength,h,dh,nbasis) 		            
C  
*------------------------------------------------------------------*
C   |  Evaluate ith basis function h and derivative dh at point z    
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)

      n = mod(i,2)
      m = i/2 + n
      a = dble(2*m-1)*2.d0*datan(1.d0)
      if (n.eq.0) then
        h = dsin(a*z)
        dh = a*dcos(a*z)/xlength
      else
        h = 1.d0 - dcos(a*z)
        dh = a*dsin(a*z)/xlength
      end if
      return
      end

C--------------------------------------------------------------------BASIS

      subroutine basis2(i,z,xlength,h,dh,nbasis)		             
C  
*------------------------------------------------------------------*
C   |  Evaluate ith basis function h and derivative dh at point z    
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)

      if (nbasis.gt.5) Stop 'Not enough basis functions in BASIS'
      go to (1,2,3,4,5),i

c.... Basis function h1(z) = z
    1 h  = z
      dh = 1.d0/xlength
      return

c.... Basis function h2(z) = 4z^2 - 3z
    2 h  = (4.d0*z - 3.d0)*z
      dh = (8.d0*z - 3.d0)/xlength
      return

c.... Basis function h3(z) = 15z^3 - 20z^2 +6z
    3 h  = ((15.d0*z - 20.d0)*z + 6.d0)*z
      dh = ((45.d0*z - 40.d0)*z + 6.d0)/xlength
      return 

c.... Basis function h4(z) = 56z^4 - 105z^3 + 60z^2 - 10z
    4 h  = ((( 56.d0*z - 105.d0)*z +  60.d0)*z - 10.d0)*z
      dh = (((224.d0*z - 315.d0)*z + 120.d0)*z - 10.d0)/xlength
      return

c.... Basis function h5(z) = z^5
    5 h  = z**5
      dh = (5.d0*z**4)/xlength
      return
     
      end
     
C-------------------------------------------------------------------APPLIED

      subroutine applied(z,dm,dq,dp,clf,n)                            
C  
*------------------------------------------------------------------*
C   |  Evaluate the distributed load functions at point z	     
 |
C   |  DM = moment, DQ = transverse force, DP = axial force	     
 |
C   |  CLF is the current load factor.				     
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      common /loads/ dmo(3),dqo(3),dpo(3)

c.... Compute the nominal values of the applied forces at point z
      f = 1.d0 - z
      go to (1,2), n

c.... Compute total transverse loads
    1 dm = dmo(1) + (dmo(2) + dmo(3)*f)*clf
      dq = dqo(1) + (dqo(2) + dqo(3)*f)*clf
      dp = dpo(1) + (dpo(2) + dpo(3)*f)*clf
      return

c.... Compute transverse loads associated with load factor only
    2 dm = dmo(2) + dmo(3)*f
      dq = dqo(2) + dqo(3)*f
      dp = dpo(2) + dpo(3)*f
      return
      end

C-------------------------------------------------------------------RESULTS

        subroutine results(x,xlength,nstep,nu,test,ndm,nbasis)
C  
*------------------------------------------------------------------*
C   |  Print results of current step to various files		     
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      logical plot
      dimension x(ndm),disp(3)
      common /out/ npts,nnstep

c.... Determine if current step is an output step
      if(mod(nstep,nnstep).eq.0) then
        plot = .true.
        nplot = npts
      else
        plot = .false.
        nplot = 1
      end if

c.... Write coefficients for current step
      if (plot) then
        write(4,2001) nstep, x(ndm), nu, test
        do i=1,nbasis
          ii = 3*(i-1)
          write(4,2002) i,(x(ii+k),k=1,3)
        end do
      end if

c.... Compute and print current geometry of beam
      if (plot) write(3,3001) nstep,x(ndm)
      z = 0.d0
      dz = 1.0/nplot
      do ii=1,nplot+1
        do j=1,3
          disp(j) = 0.d0
          do i=1,nbasis
            call basis(i,z,xlength,h,dh,nbasis)
            disp(j) = disp(j) + h*x(3*(i-1)+j)
          end do
        end do
        if (plot) write(3,3000) z*xlength,z*xlength+disp(1),disp(2)
        z = z + dz
      end do
      write(6,2003) nstep,x(ndm),(disp(j),j=1,3),nu,test
      return

 2001 format(/' Load Step  :',i4,5x,' Load factor      :',e12.5,
     *       /' Iterations :',i4,5x,' Norm of Residual :',e12.5/
     *       /'    i',15x,' a(i)',15x,' b(i)',15x,' c(i)')
 2002 format(i5,3e20.5)
 2003 format(i5,4e12.5,i5,e12.4)
 3000 format(4e15.5)
 3001 format(' Step : ',i5,5x,' Load : ',e15.5)
      end
 
C-------------------------------------------------------------------SIMPSON

      subroutine simpson(m,npoints,dz,xlength,factor)
C  
*------------------------------------------------------------------*
C   |  Evaluate the weight factor for the current integration point m
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)

      c = xlength*dz/3.d0
      n = mod(m,2)   
      if((m.eq.1).or.(m.eq.npoints)) then
         factor = c
      else if (n.eq.0) then
         factor = 4.d0*c
      else 
         factor = 2.d0*c
      endif
      return
      end

C------------------------------------------------------------------ZEROVEC

      subroutine zerovec(v,n)
C  
*------------------------------------------------------------------*
C   |  Initialize the array V to zero for all values up to n	     
 |
C  
*------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension v(n)

      do i=1,n
        v(i) = 0.d0
      end do
      return
      end
  
C--------------------------------------------------------------------INVERT

      subroutine invert(a,nmax,ndm)
C   *---------------------------------------------------------------*
C   |     Invert matrix A (nmax,nmax) when array dimension is NDM   |
C   *---------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension a(ndm,ndm)

      do 200 n = 1,nmax
        d = a(n,n)
        do 100 j = 1,nmax
          a(n,j) = -a(n,j)/d
  100   continue
        do 150 i = 1,nmax
          if(n.eq.i) go to 150
          do 140 j = 1,nmax
            if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
  140     continue
          a(i,n) = a(i,n)/d
  150   continue
      a(n,n) = 1.d0/d
  200 continue
      return
      end


