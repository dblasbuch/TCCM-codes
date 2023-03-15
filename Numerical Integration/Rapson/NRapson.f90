REAL FUNCTION NEW(X)
  REAL :: X
  NEW = ((X**2)-X)
  RETURN
END FUNCTION

MODULE DERIV

  IMPLICIT NONE

  PUBLIC :: NEW1, NEW2

  CONTAINS

  SUBROUTINE NEW1(NEW, X, h, N1) ! I redefined your first derivative as a subroutine
    ! Defining variables
    REAL :: NEW
    REAL, INTENT(IN) :: X, h
    REAL, INTENT(OUT) :: N1
    ! First derivative
    N1 = (NEW(X + h) - NEW(X)) / h
  END SUBROUTINE NEW1

  SUBROUTINE NEW2(NEW, X, h, N2) ! Apparently, passing 2 functions to a subroutine is a mess. Thus, I looked for a formula which operates everything in one step.
    ! Defining variables
    REAL :: NEW
    REAL, INTENT(IN) :: X, h
    REAL, INTENT(OUT) :: N2
    ! Second derivative
    N2 = (NEW(X + h) - 2*NEW(X) + NEW(X - h)) / h**2
  END SUBROUTINE NEW2

END MODULE DERIV

PROGRAM NR

USE DERIV, ONLY : NEW1, NEW2 ! Loading module and importing subroutines

IMPLICIT NONE

! Defining variables
REAL, EXTERNAL :: NEW
REAL :: h = 0.01
REAL :: X = 3.0
REAL :: N1, N2

WRITE(6,*) "EX2 Dani Blas Buch"
WRITE(6,*) "Calculation of the minimum of a 1D function using:"
WRITE(6,*) "The one-dimensional Newton-Raphson method"

! First check of the first derivative
CALL NEW1(NEW, X, h, N1)

DO WHILE (ABS(N1) .gt. 1E-8)

  IF (ABS(N1) .le. 1E-8 .and. ABS(X - (X - h)) .le. 1E-8) THEN

    WRITE(6,*) "The value of the minimum is:", X

  ELSE
    ! Compute both derivatives
    CALL NEW1(NEW, X, h, N1)
    CALL NEW2(NEW, X, h, N2)

    WRITE(6,*) "The value of X is:", X
    WRITE(6,*) "The value of F(X) is:", NEW(X)
    WRITE(6,*) "The value of the gradient is:", N1

    ! Newton-Raphson Method
    X = X - N1 / N2 

  ENDIF
ENDDO

END PROGRAM NR