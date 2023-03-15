     

      EXTERNAL NEW
      REAL h
      REAL X
      REAL N1, N2
      
      WRITE(6,*) "EX2 Dani Blas Buch"
      WRITE(6,*) "Calculation of the minimum of a 1D function using:"
      WRITE(6,*) "The one-dimensional Newton-Raphson method"

      h=0.01
      X=3

      CALL NEW1(NEW,X,h,N1)

      DO WHILE (ABS(N1) .gt. 1E-8)

        IF (ABS(N1) .le. 1E-8 .and. ABS(X-(X-h)) .le. 1E-8) THEN

         WRITE(6,*) "The value of the minimum is:", X

        ELSE
         
         CALL NEW1(NEW,X,h,N1)
         CALL NEW2(NEW,X,h,N2)

         WRITE(6,*) "The value of X is:", X
         WRITE(6,*) "The value of F(X) is:", NEW(X)
         WRITE(6,*) "The value of the gradient is:", N1

         X=(X-(N1/N2))

        ENDIF
      ENDDO
      END

      SUBROUTINE NEW1(NEW,X,h,N1)
      REAL NEW
      REAL X, h
      REAL N1
       
       N1=(NEW(X+h)-NEW(X))/h
      
      RETURN
      END

      SUBROUTINE NEW2(NEW,X,h,N2)
      REAL NEW
      REAL X, h
      REAL N2
       
       N2=(NEW(X+h)-2*NEW(X)+NEW(X-h))/h**2

      END
 
      FUNCTION NEW(X)
       NEW=((X**2)-X)
      RETURN
      END


