!Author: Dani Blas Buch with some help

!CODES FOR SIMPSON'S RULE and ROMBERG'S METHOD

!!    This program integrates univariable functions using simpson's rule
!!    and romberg's method. There are two subroutines accordingly, 
!!    subroutine "simpsons" and subroutine "romberg". 

!     important variables:
!     conv ==== convergence criteria (can be adjusted from main program). 
!     areas ==== integral calculated using simpson's rule, since integral
!     is a signed area under the curve. 
!     r ==== matrix containing successive integral values calculated from
!     romberg's method. 
!     msz ==== size of matrix r. 
!     ipos === a two-value array keeping the position of the converged result. 

!----------------------------------------------------------------------------------------------
! Subroutine for Simpson's rule, takes inputs such as function to be integrated (func), 
! initial number of subintervals to be employed (n), bounds of integration(a and b), 
! convergence criteria. Information that it can return via modifying values of 
! corresponding variables are: number of subintervals used (n), integral value("areas"), 
! iteration number(iter), error between two latest consecutive values(upon printing, 
! not returning in the program). When convergence is reached, loop within the subroutine stops.  
!-----------------------------------------------------------------------------------------------

      subroutine simpsons(func,n,a,b,areas,iter,conv)
      implicit double precision (a-h, o-z)         
      iter=0       ! Initial definition of iteration 
      areas_old=0  ! Arbitrary definition of area for previous iteration (for error computation). 
      error=1.0   ! Arbitrary initial definion of error to enter the loop
      write(6,*) "                                  Iterations with Simpson's rule: "
      write(6,*) " ------------------------------------------------------------------------------------------------------"          
      write(6,*) "   Iteration #     ","  Number of subintervals   ", "  Integral value  ", "                  Error  "   
      write(6,*) " ------------------------------------------------------------------------------------------------------"   
      
      do while(error.gt.conv)
      	  iter=iter+1
      	  n=n*2              ! Doubling the number of subintervals to employ (at first iteration n = 2 is used). 
          areas = 0
          h=(b-a)/(2*n)  ! Calculation of increment (x_i = a + ih, i = 0, 1..., 2n)
          
          ! Adding the contributions from endpoints         
          areas = areas + h/3 * (func(a) + func(b)) 
          
          ! Adding the contributions from odd-numbered points (i=1, 3,...,2n-1)                
          do i=1, 2*n-1, 2
              areas = areas + h/3 * func(a+i*h) * 4
          enddo
          
          ! Adding the contributions from even-numbered points.(i = 2, 4,..., 2n)    
          do i=2, 2*n-2, 2
              areas = areas + h/3 * func(a+i*h) * 2
          enddo
          
          error = dabs(areas-areas_old)      ! Calculation of error between previous and current iteration. 
          write(6,*) iter,"         ", n ,"         ", areas, "  ", error
          write(6,*) "       "
          areas_old = areas                  ! Setting new area as old area for the update for next iteration. 
      enddo
      end
  
      
!-----------------------------------------------------------------------------------
! Subroutine for Rombergs, takes inputs such as function to be integrated, bounds 
! of integration, convergence criteria. Information that it can return is value of 
! the integral at each iteration and when is the convergence criteria reached so that 
! final result can be shown ("detects" the position at which convergence is reached)
!-----------------------------------------------------------------------------------
      
      subroutine romberg(funct,a,b,conv,r,msz,ipos)
      implicit double precision(a-h,o-z)
      dimension r(1000,1000), v(1000), ipos(2)
      ! v is an array/vector keeping the values of successive h's (increments). 
      v(1)=(b-a)       ! Calculation of the very first increment     
      r(1,1) = v(1)/2*(funct(a) + funct(b)) ! Calculating first element
      do i=2,200       ! Creation of values of successive h's (increments). 
         v(i)=v(i-1)/2
      enddo
      msz = 0   ! Defining initial matrix size (can be used for gauging if there is convergence within a row). 
      i = 0  ! True/false operator for breaking the loop at the right time (see inside the loop). 
      do k = 2, 200         ! Arbitrary upper bound to calculate as much as needed to converge (200 is very high for this method). 
          term1=(r(k-1,1))  ! First term in the general formula of r(k,1)
          term2=v(k-1)      ! Second term in the general formula of r(k,1) (multiplying the third term). 
          term3=0        
          do l=1,2**(k-2)                       ! Summation to get third term of r(k,1)
              term3=term3+funct(a+(2*l-1)*v(k))
          enddo
          r(k,1)=0.5*(term1+term2*term3)        ! Combining three terms according to formula of r(k,1)
          
          do j=2,k                              ! Calculation of entire row up to diagonal (including). 
              r(k,j)=r(k,j-1) + (r(k,j-1)-r(k-1,j-1))/(4**(j-1)-1)
              
              ! If statement block to evaluate convergence each time. 
              if (dabs(r(k,j)-r(k,j-1)).lt.conv) then
                  ipos(1) = k     ! Get the first index of position of converged result 
                  ipos(2) = j     ! Get the second index of position of converged result 
                  msz = k         ! Gauge the triangular matrix size 
                  i = 1  ! If the result is converged, this will break the loop after the inner enddo statement.  
                  exit              
              endif    
          enddo
          
          if (i.ne.0) then
              do j = ipos(2)+1, k          ! Finish calculating the values of the k-th row (after which we exit the loop). 
                  r(k,j)=r(k,j-1) + (r(k,j-1)-r(k-1,j-1))/(4**(j-1)-1)
              enddo     
              exit         ! If the inner loop is broken due to convergence, outer (k dummy variable) loop will break. 
          endif              
      enddo                                            
      end


!Main program:
      
      
      implicit double precision (a-h,o-z)
      dimension r(1000,1000), ipos(2)  
      external f 
      write(6,*) " " 
      write(6,*)"This program integrates f(x)=sin(x^2)-cos(2x) from 1 to 3 with composite Simpson's rule and Romberg's method"  
      write(6,*) "------------------------------------------------------------------------------------------------------------" 
      write(6,*) "------------------------------------------------------------------------------------------------------------" 
      write(6,*) " "          
      a = 1     ! Lower bound of integration 
      b = 3     ! Upper bound of integration
      n=1       ! Initial number of points employed in simpson's rule (divided by two, will be adjusted at first iteration).  
      conv = 1D-8    ! Convergence criteria. 
      
      call simpsons(f,n,a,b,areas,iter,conv)

      write(6,*) "------------------------------------------------------------------------------------------------------" 
      write(6,*) "Final resuls with Simpson's rule after convergence: "
      write(6,*) " "
      write(6,*) "Integral value with Simpson's rule:", areas
      write(6,*) "Total number of iterations =", iter
      write(6,*) "Total number of points employed = ", 2*n + 1 
      write(6,*) "Total number of subintervals =", n   
      write(6,*) "------------------------------------------------------------------------------------------------------" 
      write(6,*) " "    ! Creating some space 
      write(6,*) " "
      write(6,*) "Results with Romberg's method: "
      write(6,*) " "
      call romberg(f,a,b,conv,r,msz,ipos)
      write(6,*) "Triangular matrix is: "
      write(6,*) " "
      
      if (msz.le.10) then ! To restrict showing the matrix up to r(10,10) 
          n=msz
      else
          n=10    
      endif
      
      do i=1,n
          write(6,200) (r(i,j), j=1,n)
      enddo       
         
      write(6,*) " "  
      write(6,*) "Converged integral result:"
      write(6,200) r(ipos(1),ipos(2))
      write(6,*) "Position of converged result in the matrix:",ipos(1),",",ipos(2) 
      write(6,*) " "
  200 format(20F18.10) ! Showing everything until 10 decimals.  

      end

           
      
!     Definition of the function to be used within the program      
      
      function f(z)
      implicit double precision(a-h,o-z)
          f = dsin(z**2) - dcos(2*z)
      end                           
      
      
      
           
      
            
