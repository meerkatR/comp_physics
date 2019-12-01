PROGRAM rootfinding
 IMPLICIT NONE
 DOUBLEPRECISION :: deltax
 DOUBLEPRECISION, PARAMETER :: deltaxfix = 5D-4, deltaxmin = 2D-5
 DOUBLEPRECISION :: x
 DOUBLEPRECISION :: xmin = -3D0, xmax = 30D0

 x = xmin
 deltax = deltaxfix

 DO WHILE (x < xmax )
!  IF ( (fcn(x)*fcn(x+deltax)) < 0D0 .AND. DABS(fcn(x)*fcn(x+deltax)) < 1D3 ) THEN
  IF ( (fcn(x)*fcn(x+deltax)) < 0D0 ) THEN
     IF (deltax < deltaxmin) THEN
      PRINT*, "root at ", x
      x = x + deltax
      deltax = deltaxfix
     ENDIF
!     else
     deltax = deltax/2D0
 !   end if
  ELSE
   x = x + deltax
  ENDIF
 END DO

 PRINT*, fcn(17.0832D0)

CONTAINS

 DOUBLEPRECISION FUNCTION fcn(x)
  IMPLICIT NONE
  DOUBLEPRECISION :: x
  
!   fcn = (x - 2D0)*(x - 17D0)*(x - 28D0)*(x - 39D0)
 ! fcn = 2D0 + 1D0/( (x-5D0) - 1D0/((x-10D0)-1d0/(x-17d0)) ) - x
 fcn = 2D0 + 1D0/( (x-5D0) - 1D0/(x-17D0) ) - x
  RETURN
  END FUNCTION

END PROGRAM


