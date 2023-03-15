!Author: Dani Blas Buch with some help


! ----------------------------------------------------------------------
! Subroutine for reading the data from file, read data is saved in
! the matrix with two columns (x,y)   
!-----------------------------------------------------------------------           

      subroutine reader(filename,p_data,length)

      implicit double precision(a-h,o-z)
      dimension p_data(200,2)
      character *20 filename
            
      open(unit=9,file=filename)
      do i = 1, 200
          read(9,*,end=25) (p_data(i,j), j = 1,2) 
          length = i       ! Once the end of the file (within rows) is reached, we gauge the length of data. 
      enddo    
   25 continue 
      close(9)
      end


! ----------------------------------------------------------------------
! Subroutine for setting up augmented matrix according to the polynomial
! regression algorithm.  
!-----------------------------------------------------------------------   
      
      subroutine set_aug_mat(p_data, length, m, aug_mat)

      implicit double precision(a-h,o-z)
      dimension p_data(200,2), aug_mat(200,200)
      
      do j=1,m+1           ! Run across rows (we need m+1 equations for fitting of polynomial of degree m). 
          el_y = 0         ! Start counting of the last column element of augmented matrix 
          do k=1,m+1       ! Run across columns (up to last column (excluding)).        
              el_x = 0     ! Starting each j,k element of the augmented matrix (excluding last column). 
              do i = 1,length   ! Length is length of the data (number of points in the data). 
                  el_x = el_x + p_data(i,1)**(j+k-2)   ! Summation of x terms in corresponding power. 
              enddo
              aug_mat(j,k) = el_x                ! Save the sum in matrix element. 
          enddo
          do i = 1,length
              el_y = el_y + p_data(i,1)**(j-1)*p_data(i,2)      ! Accumulate the last column per row . 
          enddo
          aug_mat(j,m+2) = el_y             ! Save each the last column element in row. 
      enddo               
      end           


! ----------------------------------------------------------------------
! Start of the subroutine solving systems of linear equations. 
!-----------------------------------------------------------------------   
      

      subroutine solver(a_mat_in,msize,a_mat_out)
      implicit double precision(a-h,o-z)
      dimension a_mat_in(200,200),a_mat_out(200,200),a_mat(200,200)
      nr = msize
      nc = msize + 1      
      a_mat = a_mat_in         ! Copy input matrix into the matrix for diagonalization. 
      
      write(6,*) " "
      
      ! Gaussian elimination.       
      ! Putting zeros below the main diagonal of the original matrix
      
      do i=2,nr       ! Running over the rows. 
          do j=1,i-1    ! Running over the columns (up to row number, i.e. up to diagonal). 
              cm=a_mat(i,j) ! Finding constant multiplier for the subtraction between rows. 
              do k=1,nc                 
                  a_mat(i,k) = a_mat(i,k) - cm*a_mat(j,k)   ! Subtracting the entire row. 
              enddo 
              dv = a_mat(i,j+1)                 ! Obtaining the divisor (element i,i since j runs from 1 while i runs from 2). 
          enddo
          do l=1,nc                             ! Dividing the entire row (hence, making element i,i equal to 1). 
              a_mat(i,l) = a_mat(i,l) / dv
          enddo
      enddo     
      
            
      ! ---------------------------------------------------------------
      ! Second part of the algorithm, intruducing 0's in the pivots on the left
      ! by elementary row operations. 
      ! ---------------------------------------------------------------
            
      do i=nr-1,1,-1     ! running from 1-above last row to the first row. 
          do j=nc-1,i+1,-1 ! running from last column up to the column of diagonal (not including diagonal). 
              cm=a_mat(i,j)    ! finding constant multiplier for the subtraction of rows. 
              do k=1,nc                 
                  a_mat(i,k) = a_mat(i,k) - cm * a_mat(j,k)  ! Subtracting the entire row. 
             enddo
         enddo     
      enddo                                                
      
      a_mat_out = a_mat
      end
            
      
! -----------------------------------------------------------------------------------------
! Start of the subroutine calculating R^2 parameter for determination of quality of fitting.
! R^2 = 1 - RSS/TSS (RSS === sum of squares of residuals, TSS === Total sum of squares).  
!------------------------------------------------------------------------------------------       
      
      subroutine r_squared(f,coeffs,p_data,length,m_deg,r_sq)

      implicit double precision(a-h,o-z)
      dimension coeffs(200,200), p_data(200,2)
  
      rss = 0           ! Initial definition of RSS
      do i = 1, length
          x_i = p_data(i,1)
          y_i = p_data(i,2)
          rss = rss  + (y_i - f(coeffs,m_deg,x_i))**2
      enddo                                            
      ! RSS computed. 
      
      y_sum = 0            
      do i = 1, length
          y_sum = y_sum + p_data(i,2)     ! Sum of y values computed cumulatively. 
      enddo
      y_mean = y_sum/length        ! Mean value of y values computed. 
      tss = 0
      do i = 1, length
          y_i = p_data(i,2)
          tss = tss + (y_i - y_mean)**2       ! Total sum of squares calculated according to the definition. 
      enddo
      r_sq = 1 - rss/tss               ! R^2 computed according to the definition. 
      end    
       
      
      subroutine create_output(f, coeffs,p_data,length,m_deg)

      implicit double precision(a-h,o-z)
      dimension coeffs(200,200), p_data(200,2)
      character *20 dat_plotting
      write(6,*) " "
      write(6,*) "Put the file name for data output for plotting: "
      read(5,*) dat_plotting
      
      x_ii = p_data(1,1) - 2              ! For having greater range for plotting 
      x_fi = p_data(length,1) + 2
      n_points = (x_fi - x_ii)*3000
      h_inc = (x_fi - x_ii)/n_points  
   
      open(unit=10,file=dat_plotting) 

      do i = 1, n_points + 1
          x_n = x_ii + (i-1)*h_inc
          if (i.le.length) then
              x_i = p_data(i,1)
              y_i = p_data(i,2)
              write(10,200) x_n, f(coeffs, m_deg, x_n), x_i, y_i
          else
              write(10,200) x_n, f(coeffs, m_deg, x_n)
          endif 
      enddo        
      close(10)
      write(6,*) " "
      write(6,*) "Data saved for plotting!"
  200 format(20F14.5)          
      end    
      
 
 

!Main program:    

      implicit double precision(a-h,o-z)   
      dimension p_data(200,2), aug_mat(200,200),solv_mat(200,200)
      character *20 filename, fileout
      external f
      
      write(6,*) "Provide the name of the file containing data points (x:y) in columns (no other characters):"
      read(5,*) filename 
      
      call reader(filename, p_data,length) !p_data is matrix containing points, length is size of the data. 

      write(6,*) " "
      write(6,*) "Type the polynomial degree you want to fit these data to:"
      read(5,*) m
      write(6,*) " "
      write(6,*) "           Data points "
      write(6,*) "         x             y"
      
      do i = 1,length
          write(6,200) (p_data(i,j), j = 1,2)
      enddo
      
      write(6,*) " "
      write(6,"(A42,I3)") "Degree of polynomial fitted to the data:", m
      call set_aug_mat(p_data,length,m,aug_mat)

            
      call solver(aug_mat,m+1,solv_mat)

      call r_squared(f,solv_mat,p_data,length,m, r_sq)
      
      do i = 1, m+1
          write(6,"(A38,I3,A2,F10.5) ")"Coefficient for the term with degree", i-1,":", solv_mat(i,m+2)
      enddo    

      write(6,*) " "
      write(6,"(A18,F7.5)") "R^2 value is: ", r_sq    
        
      write(6,*) " "
      write(6,*) "You can save the results of the fitting (format is degree : coefficient in two columns)"
      write(6,*) "Type name of the output file to save coefficients: "
      read(5,*) fileout

      open(unit=12,file=fileout)
      write(12,"(A18, A16)") "Degree of term", "Coefficient"
      do i = 1, m+1
          write(12,"(I12,F22.9) ") i-1, solv_mat(i,m+2)
      enddo          
      close(12)

      write(6,*) " "
      write(6,*) "Coefficients saved!"

      
      
      call create_output(f, solv_mat,p_data,length,m)      ! Creates output with degree : corresponding coefficient  
      
      write(6,*)
          
  200 format(20F14.5)          
      end
       

      
      function f(coeffs,m_deg,x)
      implicit double precision(a-h,o-z)
      dimension coeffs(200,200)
      value = 0
      do i = 1, m_deg+1
          value = value+coeffs(i,m_deg+2)*x**(i-1)
      enddo
      f = value
      end    

      
      
      
      
      
      
      
                 
      
