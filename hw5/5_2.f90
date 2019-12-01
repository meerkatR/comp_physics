!     Purpose:
!     Simulating motion of a golf ball. Include effects of: gravity (of course)
!     Calculate the ranges
!  
!     Note: we use MKS units. Meter, second, kilogram
!
!     Record of revisions:
!     Date             Programmer           Description of change
!     ===========      ==========           =====================
!     17.Oct.2009      Jen Hsu               gravity effect added.

MODULE parameters 
  IMPLICIT NONE
  SAVE
  DOUBLEPRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841971693993751
  DOUBLEPRECISION, PARAMETER :: dt = 1D-3
  DOUBLEPRECISION, PARAMETER :: mass = 0.046D0
  DOUBLEPRECISION, PARAMETER :: rho = 1.2D0, A = 0.00143D0 !air density and cross section
 
END MODULE parameters          

PROGRAM GOLF
 USE parameters
 IMPLICIT NONE

 INTEGER :: n 
 DOUBLEPRECISION :: t !the time of now
 DOUBLEPRECISION, PARAMETER :: tmax = 10D0 ! in sec
 DOUBLEPRECISION :: x, xp, vx, vxp, y, yp, vy, vyp !x for current x, xp for previous x, same for vx, y and vy
 DOUBLEPRECISION :: range 
 DOUBLEPRECISION :: theta0, v0 !initial angel, initial velocity

 OPEN(unit=14, file='5_2c.txt')

 PRINT*, "Input launching angle:"
 READ*, theta0

 x = 0D0
 y = 0D0
 v0=70d0

 vx = v0*DCOS(theta0*pi/180D0)
 vy = v0*DSIN(theta0*pi/180D0)
 
 DO n = 1, NINT(tmax/dt)
   t = n*dt 

   !x-algorithm
   xp = x
   vxp = vx
   
   yp = y
   vyp = vy

   x = xp + vxp * dt
   vx = vxp + (fgravx(xp, vxp) + fDragx(vxp, vyp) + fMagx(vxp, vyp) )/mass * dt
   !for ball without backspin
   !vx=vxp + (fgravx(xp, vxp) + fDragx(vxp, vyp) )/mass * dt


   !y-algorithm
   y = yp + vyp * dt
   vy = vyp + (fgravy(yp, vyp) + fDragy(vxp, vyp) + fMagy(vxp, vyp) )/mass * dt
   ! for ball without backspin
   !vy = vyp + (fgravy(yp, vyp) + fDragy(vxp, vyp) )/mass * dt
 
   WRITE(14,'(5F15.9)') t, x, y, vx, vy

   IF (y < 0D0) THEN
    PRINT*, "The golf ball hits the ground at t=", t , " sec,"
	PRINT*, "Range = ", x
    EXIT
   ENDIF
 
 ENDDO

 PAUSE

 CONTAINS
!External forces functions

 DOUBLEPRECISION FUNCTION fgravx(x,vx)
  IMPLICIT NONE
  DOUBLEPRECISION :: x, vx
  fgravx = 0D0
  RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION fgravy(y,vy)
  USE parameters
  IMPLICIT NONE
  DOUBLEPRECISION :: y, vy
  DOUBLEPRECISION, PARAMETER :: g = 9.8D0
  fgravy = mass*(-g)
  RETURN
 END FUNCTION
 
 DOUBLEPRECISION FUNCTION fdragx(vx,vy)
  USE parameters
  IMPLICIT NONE
  DOUBLEPRECISION :: vx, vy, v !x-component, y-component and magnitude of v
  DOUBLEPRECISION :: theta !direction of motion
  
  theta = DATAN(vy/vx)
  v = DSQRT(vx**2+vy**2)
  fDragx = -dragconst(v) * rho * A * (v**2) * DCOS(theta)
  RETURN
 END FUNCTION
 
 DOUBLEPRECISION FUNCTION fdragy(vx,vy)
  USE parameters
  IMPLICIT NONE
  DOUBLEPRECISION :: vx, vy, v 
  DOUBLEPRECISION :: theta 

  theta = DATAN(vy/vx)
  v = DSQRT(vx**2+vy**2)
  fDragy = -dragconst(v) * rho * A * (v**2) * DSIN(theta)
  RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION fMagx(vx,vy)
   USE parameters
   IMPLICIT NONE
   DOUBLEPRECISION :: vx, vy
   DOUBLEPRECISION, PARAMETER :: beta = 0.25D0 !beta = (S0 * omega/mass)
   fMagx = - beta * mass * vy
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION fMagy(vx,vy)
   USE parameters
   IMPLICIT NONE
   DOUBLEPRECISION :: vx, vy
   DOUBLEPRECISION, PARAMETER :: beta = 0.25D0 !beta = (S0 * omega/mass)
   fMagy =  beta * mass * vx
   RETURN
 END FUNCTION
  
 DOUBLEPRECISION FUNCTION dragconst(v)
   IMPLICIT NONE
   DOUBLEPRECISION, PARAMETER :: vc = 14D0 !critical transient speed
   DOUBLEPRECISION :: v
   
   !for smooth ball
   !dragconst=0.5d0
   
   IF (v < vc) THEN
     dragconst = 0.5D0
   ELSE
     dragconst = 0.5D0*vc/v
   ENDIF
 
   RETURN
 END FUNCTION
END PROGRAM
