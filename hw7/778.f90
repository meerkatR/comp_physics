PROGRAM main
  IMPLICIT NONE
  INTEGER :: seed
  DOUBLEPRECISION :: RAND
  INTEGER, PARAMETER :: N = 100
  INTEGER, PARAMETER :: maxT = 1000
  INTEGER :: step
  DOUBLEPRECISION, PARAMETER :: deltaT = 0.005D0
  INTEGER :: i,z
  DOUBLEPRECISION :: L = 10D0
  DOUBLEPRECISION, DIMENSION(N) :: xnew,xcurr,xprev !x-component of position of the i_th particle
  DOUBLEPRECISION, DIMENSION(N) :: ynew,ycurr,yprev !y-component
  DOUBLEPRECISION, DIMENSION(N) :: vx,vy !x,y-component of velocity of the i_th particle
  DOUBLEPRECISION :: v0 = 8D0
  DOUBLEPRECISION:: ax,ay !x,y-component of acceleration of the i_th particle
  DOUBLEPRECISION, DIMENSION(N,N) :: r,forcex, forcey
  DOUBLEPRECISION :: PE, KE,E  !energy
  DOUBLEPRECISION :: Vcmx,Vcmy !x,y CM velo
  DOUBLEPRECISION :: KT !Temperature
  double precision::deltaR=0.05d0
  double precision::theta0,pi


  seed = 918172
  CALL SRAND(seed)
  pi = 4D0*DATAN(1D0)

!set initial position and velocities
  Vcmx=0d0
  Vcmy=0d0
  DO i = 1, N
     xcurr(i) =  0.5d0 +dble(mod(i-1,10)) + deltar*(2D0*RAND()-1D0)
     ycurr(i) =  0.5d0 +dble(INT(i-0.01)/10) + deltar*(2D0*RAND()-1D0)
	 theta0= 2d0*pi*RAND()
	 vx(i) = v0*dcos(theta0)
	 vy(i) = v0*dsin(theta0)
	 Vcmx=Vcmx+vx(i)
	 Vcmy=Vcmy+vy(i)
  ENDDO
  Vcmx=Vcmx/dble(N)
  Vcmy=Vcmy/dble(N)

  Do i=1,N
     vx(i)=vx(i)-Vcmx
	 vy(i)=vy(i)-Vcmy
  end do


!Euler to get xprev
 DO i = 1, N
   xprev(i) = xcurr(i) - vx(i)*deltaT
   yprev(i) = ycurr(i) - vy(i)*deltaT
 ENDDO
 
DO step = 1, maxT

CALL force !calculate force before updating

DO i = 1, N
	  CALL accel(i,ax,ay)
	  xnew(i) = 2D0*xcurr(i) - xprev(i) + ax*deltaT**2
       IF ( DABS(xnew(i)) .GT. 2D0*L ) PRINT*, "TROUBLE: FOUL BALL!"
      IF ( xnew(i) .GT. L ) THEN
	    xnew(i) = xnew(i) - L
        xprev(i)= xprev(i)-L
		xcurr(i)=xcurr(i)-L
	  ELSEIF ( xnew(i) .LT. 0D0 ) THEN
	    xnew(i) = xnew(i) + L
		xprev(i)= xprev(i)+L
		xcurr(i)=xcurr(i)+L
	 ! ELSE                              
	 !   xnew(i) = xnew(i)
	!	xprev(i)= xprev(i)
      ENDIF

	  ynew(i) = 2D0*ycurr(i) - yprev(i) + ay*deltaT**2
       IF ( DABS(ynew(i)) .GT. 2D0*L ) PRINT*, "TROUBLE: FOUL BALL!"
	  IF ( ynew(i) .GT. L ) THEN
	    ynew(i) = ynew(i) - L
        yprev(i)= yprev(i)-L
		ycurr(i)=ycurr(i)-L
	  ELSEIF ( ynew(i) .LT. 0D0 ) THEN
	    ynew(i) = ynew(i) + L
		yprev(i)= yprev(i)+L
		ycurr(i)=ycurr(i)+L
	 ! ELSE                              
	 !   ynew(i) = ynew(i)
	!	yprev(i)= yprev(i)
      ENDIF
	  vx(i) = ( xnew(i) - xprev(i) )/(2d0* deltaT )
	  vy(i) = ( ynew(i) - yprev(i) )/(2d0* deltaT )

      xprev(i) = xcurr(i)
	  yprev(i) = ycurr(i)
	  xcurr(i) = xnew(i)
	  ycurr(i) = ynew(i)

    ENDDO !i

    !calculating energy, vcm, temperature
	KE = 0D0
	Vcmx = 0D0
	Vcmy = 0D0
	DO i = 1, N
      KE = KE + ( (vx(i))**2 + (vy(i))**2 )
	  vcmx = vcmx + vx(i)
	  vcmy = vcmy + vy(i)
	ENDDO
	KE=KE*0.5d0
    Vcmx=Vcmx/N
	Vcmy=Vcmy/N
	KT = KE/DBLE(N) 
	E=PE+KE
	WRITE(24,'(4ES20.10)') step*deltaT, KE,PE,E
	WRITE(25,'(4ES20.10)') step*deltaT, Vcmx, Vcmy, DSQRT(vcmx**2+vcmy**2)
	WRITE(26,'(2ES20.10)') step*deltaT, KT 

  ENDDO 

!record position
  DO i = 1, N
    WRITE(31,*) i, xnew(i), ynew(i)
	WRITE(32,*) i, vx(i)
  ENDDO

   DO i = 1, N
     DO z = 1, i-1
          write (33,*) r(i, z)
	 ENDDO
   ENDDO

 CONTAINS

 SUBROUTINE force     !calculate Fij, the force of i ON j
   IMPLICIT NONE
   INTEGER :: ifor, jfor
   DOUBLEPRECISION, PARAMETER :: rcut = 3D0
   DOUBLEPRECISION :: xtemp, ytemp  
   DOUBLEPRECISION :: costhe, sinthe  !cos and sin of theta

   PE = 0D0
   DO ifor = 1, N
     DO jfor = 1, ifor      !j <= i
       IF ( ifor == jfor ) THEN
           forcex(ifor,jfor) = 0D0
		   forcey(ifor,jfor) = 0D0
       ELSE 
         IF ( (xcurr(jfor) - xcurr(ifor)) > (L/2D0) ) THEN       !right -> left
           xtemp = xcurr(jfor) - L
         ELSEIF ( (xcurr(jfor) - xcurr(ifor)) < (-L/2D0) ) THEN  !left -> right
           xtemp = xcurr(jfor) + L
         ELSE                                !no telegraph
           xtemp = xcurr(jfor)
         ENDIF

         IF ( (ycurr(jfor) - ycurr(ifor)) > (L/2D0) ) THEN       !above -> under
           ytemp = ycurr(jfor) - L
         ELSEIF ( (ycurr(jfor) - ycurr(ifor)) < (-L/2D0) ) THEN  !under -> above
           ytemp = ycurr(jfor) + L
         ELSE                                !no telegraph
           ytemp = ycurr(jfor)
         ENDIF

         r(ifor,jfor) = distance(xcurr(ifor),ycurr(ifor),xtemp,ytemp)
		 PE = PE + LJPotential(r(ifor,jfor))
 
         IF ( r(ifor,jfor) > rcut ) THEN
           forcex(ifor,jfor) = 0D0
		   forcey(ifor,jfor) = 0D0
         ELSE
           costhe = (xtemp - xcurr(ifor))/r(ifor, jfor)   !the direction of r_ij pointing from i to j
           sinthe = (ytemp - ycurr(ifor))/r(ifor, jfor)
           forcex(ifor,jfor) = LJForce(r(ifor,jfor))*costhe
		   forcey(ifor,jfor) = LJForce(r(ifor,jfor))*sinthe
		 ENDIF 

	   ENDIF !(ifor ==jfor)
	   
     ENDDO !jfor
   ENDDO !ifor  

   DO ifor = 1, N
     DO jfor = ifor, N      !jfor >= ifor, Newton's 3rd Law
       forcex(ifor,jfor) = -forcex(jfor,ifor)
	   forcey(ifor,jfor) = -forcey(jfor,ifor)
     ENDDO !jfor
   ENDDO !ifor 
 END SUBROUTINE

 SUBROUTINE accel(i,ax,ay)  !acc of the i_th particle
   IMPLICIT NONE
   DOUBLEPRECISION, PARAMETER :: m = 1D0
   DOUBLEPRECISION :: ax, ay
   INTEGER :: i,j
   ax = 0D0
   ay = 0D0
   DO j = 1, N
     ax = ax + forcex(j,i)/m
	 ay = ay + forcey(j,i)/m
   ENDDO
   RETURN
 END SUBROUTINE

   
 DOUBLEPRECISION FUNCTION distance(xi,yi,xj,yj)
   IMPLICIT NONE
   DOUBLEPRECISION :: xi,yi,xj,yj
   distance = DSQRT( (xi-xj)**2 + (yi-yj)**2 )
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION LJForce(r)
   IMPLICIT NONE
   DOUBLEPRECISION :: r
   LJForce = 24D0 * ( 2D0 / r**13 - 1D0 / r**7)
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION LJpotential(r)
   IMPLICIT NONE
   DOUBLEPRECISION :: r
   LJpotential = 4D0*( 1D0 / r**12 - 1D0 / r**6 )
   RETURN
 END FUNCTION



END PROGRAM