!Author: Dani Blas Buch with some help
!-----------------------------------------------------------------------
!-----------NEWTON-RAPHSON METHOD---FUNCTION OPTIMIZATION---------------


! ---------------------------------------------------------------
! Start of subroutine calculating the hessian matrix (see the theory behind 
! in PDF document uploaded with source codes).     
!----------------------------------------------------------------              
      subroutine hessian(f,vec,nvar,r_hess)
      implicit double precision (a-h,o-z)
      dimension r_hess(100,100),vec(100)
      h=1D-5
      do i=1, nvar
          do j=1, nvar
              vec(i) = vec(i) + h    
              vec(j) = vec(j) + h
              term1 = f(vec)                   ! vec = (x1,...,x_i+h,...,x_j + h,...,x_n)
              vec(i) = vec(i) - 2*h
              term2 = f(vec)                   ! vec = (x1,...,x_i-h,...,x_j + h,...,x_n)
              vec(j) = vec(j) - 2*h
              term4 = f(vec)                    ! vec = (x1,...,x_1-h,...,x_j - h,...,x_n)
              vec(i) = vec(i) + 2*h
              term3 = f(vec)                   ! vec = (x1,...,x_i+h,...,x_j - h,...,x_n)
             
              ! Calculation and saving i,j element of the hessian.              
              r_hess(i,j) = 0.25*(term1-term2-term3+term4)/h**2
              
              ! Recovery of the vector for next iteration 
              vec(i) = vec(i) - h
              vec(j) = vec(j) + h  ! vec = (x1,...,x_i,x_j,...x_n)              
          enddo 
      enddo        
      end             
  
      
      
! ---------------------------------------------------------------
! Start of the subroutine for calculating the gradient vector. 
!----------------------------------------------------------------                      
      subroutine gradient(f,vec,nvar,grad)
      implicit double precision (a-h, o-z)
      dimension grad(100), vec(100)
      h=1D-6
      do i = 1, nvar
          term2 = f(vec)       ! vec = (x1,...,x_i,..., x_n)
          vec(i) = vec(i) + h  ! vec = (x1,...,x_i+h,..., x_n)
          term1 = f(vec)      
          ! Calculation of the i-th element of the gradient. 
          grad(i) = (term1 - term2)/h 
          ! Recovery of the vector for next iteration 
          vec(i) = vec(i) - h
      enddo
      end 
      
   
      
! ---------------------------------------------------------------
! Start of the exchange subroutine to use within inversion subroutine.
!----------------------------------------------------------------                
      subroutine exchange(n,m,p)
      implicit double precision (a-h,o-z)
      dimension p(100,100), b(100,100), s(100)
          
      ! p is the matrix submitted to be treated for exchanges.            
      ! m is the row/column index at which 0 was found (diagonal element). 
      ! n is size of the matrix. 
          
      b = p        ! b is another matrix for manipulation (auxiliary). 
      
      ! This goes through every row below the one that has 0 in diagonal 
      ! and as soon as finds the row which is appropriate to exchange with,
      ! it does so (and breaks the loop immediately). 
      do i = m+1,n 
          if (b(i,m).ne.0) then
              s = b(m,:)    ! s is auxiliary variable to keep entire row
              b(m,:) = b(i,:)
              b(i,:) = s
              p=b          ! Copy the matrix back to p. 
              exit
          endif
      enddo                
      end      


! Start of the matrix inversion subroutine.    
            
      subroutine a_mat_inv(nvar,or_mat,ainv_mat)     
      implicit double precision (a-h,o-z)
      dimension or_mat(100,100), ainv_mat(100,100), p(100,100)
      p = or_mat   ! Creating the copy of original matrix (input matrix). 
     
      ! nvar is number of variables in the function, defining size of
      ! the hessian matrix. ainv_mat is inverted matrix. 
      ! Create partitioned n x 2n matrix for inversion algorithm. 
      ! (Right hand side being the identity matrix, left - the same). 
      
      do i=1,nvar
          do j=n+1,2*nvar
              if (i+nvar.eq.j) then
                  p(i,j) = 1
              else
                  p(i,j) = 0
              endif
          enddo
      enddo     

      
      ! --------------------------------------------------------------
      ! Treating the matrix if it has any diagonal elements 0. 
      do i = 1,n 
          if (p(i,i).eq.0) then
              call exchange(nvar,i,p)
          endif 
      enddo    
      
           
      ! --------------------------------------------------------------
      ! Start of the Gaussian elimination procedure (reducing to upper
      ! triangular matrix form the left part of the nvar x 2nvar matrix).  
      ! -------------------------------------------------------------
      
      ! Making p(1,1) = 1 and adjusting the first row accordingly
      dv=p(1,1)        
      do j=1,2*nvar
          p(1,j) = p(1,j)/dv 
      enddo            

      ! Gaussian elimination.       
      ! Putting zeros below the main diagonal of the original matrix
      
      do i=2,nvar       ! Running over the rows. 
          do j=1,i-1    ! Running over the columns (up to row number, i.e. up to diagonal). 
              cm=p(i,j) ! Finding constant multiplier for the subtraction between rows. 
              do k=1,2*nvar                 
                  p(i,k) = p(i,k) - cm*p(j,k)   ! Subtracting the entire row. 
              enddo 
              dv = p(i,j+1)
          enddo
          do l=1,2*nvar
              p(i,l) = p(i,l) / dv
          enddo
      enddo
      
      
      ! ---------------------------------------------------------------
      ! Second part of the algorithm, intruducing 0's in the pivots on the left
      ! by elementary row operations. 
      ! ---------------------------------------------------------------
            
      do i=nvar-1,1,-1     ! running from 1-above last row to the first row. 
          do j=nvar,i+1,-1 ! running from last column up to the column of diagonal (not including diagonal). 
              cm=p(i,j)    ! finding constant multiplier for the subtraction of rows. 
              do k=1,2*nvar                  
                  p(i,k) = p(i,k) - cm * p(j,k)  ! Subtracting the entire row. 
             enddo
         enddo     
      enddo                                     
      
      
      !-------------------------------------------------------------------
      ! Taking right-hand part of the matrix, and saving as and inverse matrix.
      do i=1,nvar  ! running accross rows 
          do j = nvar+1,2*nvar ! running across columns only to the right side of p(i,j). 
              ainv_mat(i,j-nvar) = p(i,j)   ! saving right hand side in the new matrix. 
          enddo
      enddo
 

! --------------------------------------------------------------
! SUBROUTINE FOR MULTIPLYING MATRIX BY THE VECTOR
! --------------------------------------------------------------

      subroutine p_mat_vec(a_mat,nvar,vect_in, vect_out)
      implicit double precision (a-h,o-z)
      dimension a_mat(100,100), vect_in(100), vect_out(100)
      
      do i=1,nvar          ! Run over rows
          el=0             ! Element counting started
          do j=1,nvar      ! Run over columns of matrix
              el = el + a_mat(i,j)*vect_in(j) ! Multiply element by element
          enddo
          vect_out(i) = el     ! Save obtained element in new vector element
      enddo      
      end

!Main program: 

      implicit double precision (a-h,o-z)
      dimension r_hess(100,100),x(100),grad(100),ainv_mat(100,100),x_new(100), h_gr(100),x_old(100)
      external f    
      write(6,*) " "            ! Create initial space in output 
      write(6,*) "----------------------------------------------------------------------------"  
      write(6,*) "----------------------------------------------------------------------------"        
      write(6,*) "This program optimizes functions using newton-raphson method" 
      write(6,*) "In this case, we are optimizing f(x) = sin(x+y) + (x-y)^2 - 1.5x + 3.5y + 3"      

      nvar = 2           ! Number of variables in our function. 
      x(1) = 1.0         ! Starting point of coordinate 1
      x(2) = 3.0         ! Starting point of coordinate 2
      tol_error = 1D-8   ! Error tolerance for the gradient and distance between two points. 
      x_new(1) = 2.0     ! Arbitrary definition
      x_new(2) = 4.0     ! Arbitrary definition. 
      ! Tolerance will be considered against norm of the gradient and
      ! against the norm of the vector difference between two points (the distance)
      x_old(1) = 0.1
      x_old(2) = 0.1
      distance = 1        ! arbitrary definition to enter the loop (distance between successive points)    
      iter = 1            ! Iteration counter  
      error = 1           ! Arbitrary definition of one of the criterion to enter the loop (error between two successive values). 
      write(6,*) "----------------------------------------------------------------------------"  
      write(6,*) "----------------------------------------------------------------------------"    
      
      do while ((error.gt.tol_error).and.(distance.gt.tol_error))
          write(6,*) " "
          write(6,*) "Iteration: #",iter     ! Print iteration #. 
          write(6,*) " "                     ! Space for better seeing the results. 
          
          ! Writing values of variables and value of function at that point.   
          write(6,*) "      x           y" 
          write(6,200) x(1), x(2)                      !x variable is x(1) and y variable is x(2)
          write(6,*) " "
          write(6,"(A12,F12.6)") "     f(x,y):", f(x,nvar) 
          write(6,*) " "
                 
          call gradient(f,x,nvar,grad)         ! computing gradient
          error = dabs(f(x,nvar) - f(x_old,nvar))    ! Calculation of error based on old and current values of function. 
          ! Writing gradient on the screen
          write(6,"(A14,100F12.8)")"Gradient:", (grad(i),i=1,nvar)  
          write(6,*) " "        
          
          !computing hessian matrix          
          call hessian(f,x,nvar,r_hess)   
          
          ! Writing hessian on the screeen.           
          write(6,*) "   Hessian:"
          do i=1,nvar
              write(6,200) (r_hess(i,j), j=1,nvar)
          enddo                
          
          call a_mat_inv(nvar,r_hess,ainv_mat) !computing inverse of hessian                     
          call p_mat_vec(ainv_mat,nvar,grad,h_gr) ! computing product of inverse of hessian and gradient
          
          ! h_gr ==== product of inverse of hessian (ainv_mat) and gradient(grad). 
          ! ------------------------------------------------------------
          ! Start of calculation of the new point     
          
          do i=1,nvar
              x_new(i) = x(i) - h_gr(i)
          enddo   
           
          
          !-------------------------------------------------------------
          !Start of calculation of distance between old and new points
                    
          distance = 0
          do i = 1,nvar
              distance = distance + (x_new(i)-x(i))**2
          enddo
          distance = dsqrt(distance)   ! Taking square root to calculater true distance.
          
          ! End of calculation of distance between old and new points---
          !-------------------------------------------------------------
          x_old = x
          x = x_new     ! Updating the current point
          write(6,*) " "
          write(6,"(A,F12.10)") " Convergence in distance between points: ", distance
          write(6,"(A,F12.10)") " Convergence in value of function: ", error 
          write(6,*) " "  
          write(6,*) "-----------------------------------------------------------------"          
          iter = iter+1    ! Increase the iteration by 1. 
      enddo            
      !---------------------------------------------------------------
      ! End of the while loop
      ! --------------------------------------------------------------
      write(6,*) "----------------convergence criteria satisfied-------------------"
      write(6,*) "-----------------------------------------------------------------"    
      write(6,*) " "
      write(6,*) " "
      write(6,*) "Final results:"
      write(6,*) " "
      write(6,*) "Number of iterations needed:", iter-1 !(Because just before exiting while loop we increase it by 1)
      write(6,*) " "
      write(6,*) "Coordinates of the minimum:"
      write(6,*) "      x           y" 
      write(6,200) x(1), x(2)                      !x variable is x(1) and y variable is x(2)    
      write(6,*) " "  
      write(6,"(A14,100F12.8)")"Gradient:", (grad(i),i=1,nvar) 
      write(6,*) " "              
      write(6,*) "   Hessian:"
      do i=1,nvar
          write(6,200) (r_hess(i,j), j=1,nvar)
      enddo
      write(6,*) " "                         
  200 format(20F12.6)              
      end 
     
      function f(x,nvar)
      implicit double precision(a-h,o-z)
      dimension x(nvar)
      f = dsin(x(1)+x(2)) + (x(1)-x(2))**2 - 1.5*x(1) + 3.5*x(2) + 3
      end            


