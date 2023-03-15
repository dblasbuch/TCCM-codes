!Author: Dani Blas Buch with some help
!FUNCTION INTEGRATION (Gauss Legendre Method)

!!    Main variables:
!!    n === number of quadrature points. 
!!    t === vector containing Gauss-Legendre zeros (points to employ) for particular n
!!    w === array/vector containing weights for corresponding points
!!    X1,X2 === bounds of integration. 
!!    XM,XL === parameters (computed) for variables change. 
!!    par_m and par_c === corresponding parameters for change of variables.


      SUBROUTINE Sub_GauLeg(X1,X2,t,w,n)

      IMPLICIT NONE
      
      INTEGER :: m,i,j
      INTEGER,intent(in) :: n   !Number of Gaussian points
      REAL*8, dimension(n),intent(out) :: W,t
      REAL*8 :: XM,XL,X1,X2,EPS,P1,P2,P3,pi,Z1,Z,PP

	!Relative precision
   	EPS = 1.D-8

  	!double precision arccosine. Pi=3.14159
 	pi = DACOS(-1.D0)

	!N = number of Gauss Points
	!Roots are symmetric in the interval (find half of them)  
   	m = (n + 1) / 2
	
	!Coats= -1 and 1, Gauss-Legendre 
	!Variable change
      	XM=0.5D0*(X1+X2)
      	XL=0.5D0*(X2-X1)


	!Loop over the desired roots
      	
         DO i = 1,m
         Z = DCOS (pi * (i - 0.25D0)/(n + 0.5D0))

	!Starting with the above approximation to the i-th root,
	!loop of refinement by NEWTON'S method   
10      	P1 = 1.D0
        	P2 = 0.D0

		!Loop to get the Legendre polynomial evaluated at z                
         	DO j = 1,n
           	  P3 = P2
           	  P2 = P1
            	  P1 = ((2.D0 * j - 1.D0) * Z * P2 - (j - 1.D0) * P3)/j
         	END DO
!p1 = the Legendre polynomial we want.
!Compute pp, its derivative, using the polynomial of one lower order.
         PP = n * (Z * P1 - P2)/(Z * Z - 1.D0)
         Z1 = Z
         Z = Z1 - P1/PP	      ! Newton's Method  */

         IF (DABS(Z-Z1) .GT. EPS) GO TO 10

	! Roots (symmetric) 
         t(i) = XM - XL * Z
         t(n + 1 - i) = XM + XL * Z
	!Compute the weight and its symmetric counterpart 
         W(i) = 2.D0 * XL/((1.D0 - Z * Z) * PP * PP)
         W(n + 1 - i) = W(i)

         END DO  

         END SUBROUTINE Sub_GauLeg      
      
      
      

!Main program:  
      
      Implicit double precision (a-h,o-z)

      dimension t(1000), w(1000)   ! Arbitrary dimension definitions. 
      external f 
      write(6,*) " "
      
      write(6,*)"This program integrates f(x)=sin(x^2)-cos(2x) using Gauss-Legendre method"  !Info to user
      write(6,*)" "

      a = 1
      b = 3      

      write(6,"(A29, F7.4)") "Lower bound of integration:", a       
      write(6,"(A29, F7.4)") "Upper bound of integration:", b

      write(6,*)" "
      write(6,*)" "

      write(6,*) "Iterations with Gauss-Legendre method:"

      write(6,*) " "
      write(6,*) " "

      conv = 1D-8      ! Convergence criterion selected. 
      g_error = 1.0    ! Arbitrary error definition to enter while loop. 
      par_m = (b-a)/2  ! Parameter m in the calculation. (See pdf file for explanation). 
      g_area_old = 0   ! Arbitrary definition of old area. 
      n=2              ! Initial number of quadrature points 
      iter = 0         ! Iteration defined to be increased at each step. 
     
      do while((g_error.gt.conv).and.(n.lt.11))
      g_area=0

      call Sub_GauLeg(a,b,t,w,n)    
          iter = iter + 1
          do i = 1, n
          g_area = g_area + w(i)*f(t(i))
      enddo

          g_area = g_area * par_m          !multiply by factor to match the integral. 
          g_error = dabs(g_area_old-g_area)
          g_area_old = g_area
          write(6,*) "Iteration:", iter
          write(6,*) "Number of quadrature points:", n
          write(6,*) "Integral value:", g_area
          write(6,*) "Absolute Error:", g_error
          write(6,*) " "
          n = n+1   
      enddo

      write(6,*) " "
      write(6,*) "Final results:"
      write(6,*) "The integral value is:", g_area, "Number of quadrature points:", n-1  
      write(6,*) " "                   
      end  
      
      
!Definition of the function. 

      function f(z)
      implicit double precision (a-h,o-z)
           f = dsin(z**2) - dcos(2*z)
      end      
      
      

