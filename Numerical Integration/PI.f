      INTEGER S

      WRITE(6,*) "EX2 Dani Blas Buch"
      WRITE(6,*) "Calculation of pi number"

      H=10
      CIRC=0
      OT=0
      TOT=0
      IE=10
      
      L=RAN(M)
      J=RAN(K)

      DO S=0,IE,10
       DO I=1+S,10+S
        X=RAN(L)
        Y=RAN(J)
         K=1+(H*X)
         N=1+(H*Y)

         WRITE(6,*) "The random numbers are:", K,"and", N

        DIST=SQRT((((H/2)-K)**2)+(((H/2)-N)**2))

        IF (DIST .lt. H/2) THEN
         CIRC=CIRC+1
          TOT=TOT+1
        ELSE
         OT=OT+1
         TOT=TOT+1
        ENDIF

        PI=((4*CIRC)/TOT)

        ERR=(3.1415926535898-PI)

        ENDDO

       WRITE(6,*) "The total number of points is:", TOT
       WRITE(6,*) "The number of points inside the circle is:", CIRC
       WRITE(6,*) "The calculated pi number is:", PI
       WRITE(6,*) "The error in comparison to the real number is:", ERR
      
       ENDDO
  
       END

       
      
       
