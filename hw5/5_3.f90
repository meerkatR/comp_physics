!     Purpose:
!     Simulating motion of a pendulum. With RK4 algorithm
!  
!     Note: we use MKS units. Meter, second, kilogram
!
!     Record of revisions:
!     Date             Programmer           Description of change
!     ===========      ==========           =====================
!     12.Oct.2009      Jen Hsu   

MODULE parameters
  IMPLICIT NONE
  SAVE
  DOUBLEPRECISION, PARAMETER :: pi = 3.14159265358979323846264338328
  DOUBLEPRECISION, PARAMETER :: fd = 1.2D0 ! co-eff of driving force term
  DOUBLEPRECISION, PARAMETER :: omegaD =2d0/3d0    !0.66666666666666666666D0
  DOUBLEPRECISION, PARAMETER :: tp = 3d0*pi !period=(2 Pi/omegaD)
  DOUBLEPRECISION, PARAMETER :: dt = tp/300d0 !so we can keep track at time integer multiple of tp
 END MODULE parameters           

PROGRAM pendulum
 USE parameters
 IMPLICIT NONE
 !INTEGER :: ierr
 INTEGER :: i
 DOUBLEPRECISION :: t
 DOUBLEPRECISION, PARAMETER :: tmax = 5000D0 ! sec
 DOUBLEPRECISION :: theta, omega
 DOUBLEPRECISION, DIMENSION(INT(tmax/dt)) :: thetaarray
 DOUBLEPRECISION :: frequency

 OPEN(unit=14, file='data5_3a.txt')
 OPEN(unit=21, file='data5_3b.txt')
 OPEN(unit=15, file='data5_3c.txt')


 theta = 0.2D0
 omega = 0D0
 
 DO i = 1, INT(tmax/dt)
   t = i*dt
   CALL rk4(t, theta, omega, fx, fvx)
  
   !make theta periodic
   IF (theta > pi) THEN
     theta = theta - 2D0*pi
   ELSEIF (theta < -pi) THEN
     theta = theta + 2D0*pi
   ELSE
     theta = theta
   ENDIF
   thetaarray(i) = theta   

   WRITE(14,'(5F15.9)') t, theta, omega

!data for Poincare map
   IF ( MOD(i, 300) == 0) THEN  !record the point at t=n*period
    WRITE(21, *) t, theta, omega
   ENDIF
 ENDDO


!print data for power spectrum

 DO frequency = 0D0, 1D0, 0.001D0
   WRITE(15, *) frequency, powerSpect(thetaarray, frequency)
 ENDDO

 PAUSE

 CONTAINS

!the differentiation equation for theta d theta/dt = omega
 DOUBLEPRECISION FUNCTION fx(t, theta, omega)
   USE parameters
   IMPLICIT NONE
   DOUBLEPRECISION :: t, theta, omega
   fx = omega
   RETURN
 END FUNCTION
   

!the differentiation equation for omega d omega/dt = fvx
 DOUBLEPRECISION FUNCTION fvx(t, theta, omega)
   USE parameters
   IMPLICIT NONE
   DOUBLEPRECISION :: t, theta, omega
   DOUBLEPRECISION :: fHook, fdrag, fdrive

   fvx = -DSIN(theta)-0.5d0*omega+fd* DSIN(omegaD * t)
   RETURN
 END FUNCTION

 SUBROUTINE rk4(t, x, vx, fx, fvx)
   USE parameters
   EXTERNAL fx, fvx
   DOUBLEPRECISION :: fx, fvx
   DOUBLEPRECISION :: t, x, vx  ! x and vx are both inputs and outputs. the values will be updated
   DOUBLEPRECISION :: k1x, k2x, k3x, k4x
   DOUBLEPRECISION :: k1vx, k2vx, k3vx, k4vx

   k1x  = fx(t, x, vx) * dt
   k1vx = fvx(t, x, vx)* dt

   k2x  = fx(t+0.5D0*dt, x+0.5D0*k1x, vx+0.5D0*k1vx) * dt
   k2vx = fvx(t+0.5D0*dt, x+0.5D0*k1x, vx+0.5D0*k1vx)* dt

   k3x  = fx(t+0.5D0*dt, x+0.5D0*k2x, vx+0.5D0*k2vx) * dt
   k3vx = fvx(t+0.5D0*dt, x+0.5D0*k2x, vx+0.5D0*k2vx)* dt

   k4x  = fx(t+dt, x+ k3x, vx+k3vx) * dt
   k4vx = fvx(t+dt, x+ k3x, vx+k3vx) * dt

   x = x + (k1x + 2D0*k2x  + 2D0*k3x  + k4x)/6D0
   vx=vx + (k1vx+ 2D0*k2vx + 2D0*k3vx + k4vx)/6D0
     
 END SUBROUTINE
 
 DOUBLEPRECISION FUNCTION powerSpect(array,f)
   USE parameters
   IMPLICIT NONE
   INTEGER :: i
   DOUBLEPRECISION :: f
   DOUBLEPRECISION :: FReal, FIm
   DOUBLEPRECISION, DIMENSION(INT(tmax/dt)) :: weight
   DOUBLEPRECISION, DIMENSION(INT(tmax/dt)) :: array

!Integration by trapzoid method
  !generating the weights
  DO i = 1, INT(tmax/dt) !dtermining tpart and weight
     IF (i ==1 .OR. i == INT(tmax/dt)) THEN
       weight(i) = 0.5D0
     ELSE
       weight(i) = 1D0
     ENDIF
  ENDDO

 !Real part integral
  FReal = 0D0
  DO i = 1, INT(tmax/dt)
   FReal = FReal + array(i)*DCOS(pi*f*(i*dt))*weight(i)
  ENDDO 
   FReal = FReal*dt
   
  !Imaginary part integral
  FIm = 0D0
  DO i = 1, INT(tmax/dt)
   FIm = FIm + array(i)*DSIN(pi*f*(i*dt))*weight(i)
  ENDDO !tindex
   FIm = FIm*dt
   powerSpect = FReal**2 + FIm**2
  RETURN
END FUNCTION

END PROGRAM

