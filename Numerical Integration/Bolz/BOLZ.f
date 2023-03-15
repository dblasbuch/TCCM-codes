

      EXTERNAL DIM

      WRITE(6,*) "EX1 Dani Blas Buch"
      WRITE(6,*) "Calculation of roots by Bolzano's methods:"
      WRITE(6,*) "a) The Bisection method"
      WRITE(6,*) "b) The Regula Falsi method"

      a=3
      b=4

      IF (DIM(a) .lt. 0 .and. DIM(b) .gt. 0) THEN
       
       CALL BIS(DIM,a,b,c,rebis)
       CALL FAL(DIM,a,b,d,refal)

       WRITE(6,*) "The roots using the Bisection method for are:", rebis
       WRITE(6,*) "The roots using the Regula Falsi method for are:", refal
      
      ELSE 

       WRITE(6,*) "Introduce two other values for the calculation"
       WRITE(6,*) "The previous ones were", a, b
       WRITE(6,*) "a="
       READ(5,*) a
       WRITE(6,*) "b="
       READ(5,*) b

       CALL BIS(DIM,a,b,c,rebis)
       CALL FAL(DIM,a,b,d,refal)

       WRITE(6,*) "The roots using the Bisection method are:", rebis
       WRITE(6,*) "The roots using the Regula Falsi method are:", refal

      ENDIF

      END



      SUBROUTINE BIS(DIM,a,b,c,rebis)

      c=((a+b)/2)

      DO WHILE (DIM(c) .gt. 1E-8)

      c=((a+b)/2)
 
      IF (DIM(c) .gt. 1E-8) THEN
       WRITE(6,*) "The root is in the interval (a,c)"
       I=c
       b=I

      ELSE

       WRITE(6,*) "The root is in the interval (c,b)"
       J=c
       a=J

      ENDIF

      ENDDO

      rebis=c

      END



      SUBROUTINE FAL(DIM,a,b,d,refal)

       d=((a*DIM(a))-(b*DIM(b)))/(DIM(b)-DIM(a))

      DO WHILE (DIM(c) .gt. 1E-8)

       d=((a*DIM(a))-(b*DIM(b)))/(DIM(b)-DIM(a))
      
      IF (DIM(d) .gt. 1E-8) THEN
       WRITE(6,*) "The root is in the interval (a,d)"
       H=d
       f=H

      ELSE

       WRITE(6,*) "The root is in the interval (d,b)"
       K=d
       e=K

      ENDIF

      ENDDO

      refal=d

      END
 

       
      FUNCTION DIM(X)
       e=2.718281828459
       DIM=(e**(-X))*(3.2*(SIN(X)-((1/2)*COS(X))))
      RETURN
      END
      
