!     Purpose:
!     Bifurcation diagram
!  
!    
!
!     Record of revisions:
!     Date             Programmer           Description of change
!     ===========      ==========           =====================
!     15.Oct.2009      Jen Hsu   

PROGRAM logisticsmap
 IMPLICIT NONE

 INTEGER :: i

 DOUBLEPRECISION :: x
 DOUBLEPRECISION :: mu !growth parameter

 OPEN(unit=14, file='data5_4.txt')

 DO mu = 0D0, 4D0, 0.005D0
! for mu=3.8 and 4.0 case
!   print *, "Input mu:"
!   read(*,*) mu
   x = 0.6D0
   DO i = 1, 1000  !Wait until the dust settles
     x = func(mu, x)
   ENDDO !xiter    

   DO i = 1, 2000 !start to print the value of mu and x
     x = func(mu, x)
	 WRITE(14, '(F8.5)') x
   ENDDO 
   print *, "mu=", mu    

 ENDDO

 PAUSE

 CONTAINS
 
 DOUBLEPRECISION FUNCTION func(mu, x)
  IMPLICIT NONE
  DOUBLEPRECISION :: x, mu
  func = mu * x * (1D0 - x)
  RETURN
 END FUNCTION


END PROGRAM
