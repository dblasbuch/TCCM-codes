!Author: Dani Blas Buch with some help
!-----------STEEPEST DESCENT METHOD---FUNCTION OPTIMIZATION--------------
     
!!    Important variables:
      
!!    nvar ==== number of variables in the function, dimension of gradient vector. 
!!    distance ==== distance between two successive points upon optimization
!!    (calculated as square root of sum over (x_i(new) - x_i("old"))^2 + ... + (x_n(new) - x_n("old"))^2 

!!    vnorm ==== norm (length) of the gradient vector


! Subroutine for calculating the gradient vector: 
                    
      subroutine gradient(f,vec,nvar,grad)
      implicit double precision (a-h, o-z)
      dimension grad(nvar), vec(100)
      h=1D-6
      do i = 1, nvar
          term2 = f(vec)       ! vec = (x1,...,x_i,..., x_n)
          vec(i) = vec(i) + h  ! vec = (x1,...,x_i+h,..., x_n)
          term1 = f(vec)              
          grad(i) = (term1 - term2)/h    ! Calculation of the i-th element of the gradient (finite difference method). 
          vec(i) = vec(i) - h    ! Recovery of the vector for next iteration, ! vec = (x1,...,x_i,..., x_n) 
      enddo
            
      ! Normalization of gradient
      vnorm = dsqrt(dot_product(grad,grad))   ! Calculation of the norm of gradient. 
      do i=1, nvar
          grad(i) = grad(i)/vnorm      ! Normalization of the gradient      
      enddo
      
      end      
      
!Main program:  
      
      implicit double precision(a-h,o-z)
      dimension coord(100,100), grad(100)
      external f
      write(6,*) " "
      write(6,*) "----------------------------------------------------------------------------"  
      write(6,*) "----------------------------------------------------------------------------"         
      write(6,*) "This program optimizes functions using steepest (gradient) descent method" 
      write(6,*) "In this case, we are optimizing f(x) = sin(x+y) + (x-y)^2 - 1.5x + 3.5y + 3"  
      iter = 1      ! initial iteration number 
      gm = 0.3      ! Step size which is fixed 
      nvar = 2   ! Number of variables within the function -- Which can be adjusted according to which function we want to optimize. 
      distance = 1 ! (In principle this would be the convergence criteria: P_(k+1)-P_k ≤ tolerance). 
                     ! (However, this at all iterations (except first) will be equal to the step size, which is constant for us).       
      coord(1,1) = 1   
      coord(1,2) = 3  ! Starting points for the minimization (x and y, respectively). 
      write(6,*) "----------------------------------------------------------------------------"  
      write(6,*) "----------------------------------------------------------------------------"  
                 
      
      !Steepest Descent Algorithm:

      do while (iter.lt.41)
          write(6,*) " "      
          call gradient(f,coord(iter,:),nvar,grad) 
          
          ! coord(iter,:) is a vector of nvar dimensions representing point for each iterations.           
          write(6,*) "Iteration: #",iter
          write(6,*) " "
          write(6,*) "      x           y"  
          write(6,198) coord(iter,1),coord(iter,2)              ! Write values of x and y on the screen 
          write(6,*) " "
          write(6,"(A,F12.6)")  "    f(x,y): ", f(coord(iter,:),nvar)   
          write(6,*) " "
          write(6,"(A,A1,F10.6,A2, F10.6,A2)") "   Gradient: ", "(",grad(1),", ", grad(2)," )"
          write(6,*) " "
          distance = 0
          do i = 1,nvar
              distance = distance + (coord(iter,i)-coord(iter-1,i))**2    ! Calculation of distance^2 between successive points. 
          enddo    
          distance = dsqrt(distance)       ! Taking square root to calculate the distance (note that one variable name is used). 
           
          val_old = f(coord(iter-1,:),nvar)     ! Calculation of function value at previous point
          val_new = f(coord(iter,:),nvar)       ! Calculation of function value at current point .
          error=dabs(val_new - val_old)         ! Calculation of absolute difference between successive values of function.            
  
          write(6,"(A,F12.6)") "Error in value of function: ", error
          write(6,"(A,F12.6)") "Error in distance between points:", distance
          write(6,*) " "         
          
          iter = iter + 1     ! Increase the iteration number 
          

          do i = 1,nvar                     
              coord(iter,i) = coord(iter-1,i) - gm*grad(i)   ! Calculation of the next point 
          enddo    
          
  
      enddo
         
      write(6,*) "Final results:"
      write(6,*) "Total number of iterations:", iter-1 !(Because just before exiting while loop we increase it by 1)
      write(6,*) " "
      write(6,*) "Coordinates at final iteration:"
      write(6,*) "      x           y" 
      write(6,198) coord(iter-1,1),coord(iter-1,2)
      write(6,*) " "
      write(6,"(A, F12.6)") "   f(x,y) at final iteration:  ", f(coord(iter-1,:),nvar)
      write(6,*) " " 
      
      write(6,"(A,A1,F10.6,A2,F10.6,A2)") "   Gradient at final iteration: ", "(",grad(1),", ", grad(2)," )"
      write(6,*) " "
      write(6,*) "-----------------------------------------------------------------"              
      write(6,*) "-----------------------------------------------------------------"  
    !Representing successive (x,y) coordinates upon approaching minimum. 
      write(6,*) "Successive coordinates upon approaching minimum are:"
      write(6,*) "-----------------------------------------------------------------"  
      write(6,*) "       x           y"	  
      do i=1,iter-1
          write(6,198) (coord(i,j), j=1,2) 
      enddo
  198 format(10F12.6)     ! For rounding the values where convenient. .              
      end
       
   
      function f(p,nvar)
      implicit double precision(a-h,o-z)
      dimension p(nvar)
      f = dsin(p(1)+p(2)) + (p(1)-p(2))**2 - 1.5*p(1) + 3.5*p(2) + 3
      end            




      

