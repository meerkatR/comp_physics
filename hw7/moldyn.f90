PROGRAM moleclar
!  Purpose:
!   molecular dynamics
!   Verlet method
!
!
!  Record of revisions:
!     Date             Programmer           Description of change
!     ====             ==========           =====================
!   18.Nov.2009        Jen                  Kick-off
!

  IMPLICIT NONE
  
  INTEGER :: seed
  DOUBLEPRECISION :: RAND

  INTEGER, PARAMETER :: N = 100
  INTEGER, PARAMETER :: maxT = 5000

  INTEGER :: it

  DOUBLEPRECISION, PARAMETER :: deltaT = 0.005D0
  
  INTEGER :: i,j
  
  DOUBLEPRECISION :: sidelen = 10D0
  
  DOUBLEPRECISION, DIMENSION(N) :: px !x-component of position of the i_th particle
  DOUBLEPRECISION, DIMENSION(N) :: py !y-component
  DOUBLEPRECISION, DIMENSION(N) :: px_old, px_older !
  DOUBLEPRECISION, DIMENSION(N) :: py_old, py_older !
  
  DOUBLEPRECISION, DIMENSION(N) :: vx !x-component of velocity of the i_th particle
  DOUBLEPRECISION, DIMENSION(N) :: vy !y-comp
  DOUBLEPRECISION :: v0 = 0D0

  DOUBLEPRECISION, DIMENSION(N) :: ax !x-component of acceleration of the i_th particle
  DOUBLEPRECISION, DIMENSION(N) :: ay !y-component
  
  DOUBLEPRECISION, DIMENSION(N,N) :: forcex, forcey
  DOUBLEPRECISION :: pe, ke  !potential and kinetic energy
  DOUBLEPRECISION :: vcmx !x CM velo
  DoUBLEPRECISION :: vcmy !y CM velo
  DOUBLEPRECISION :: temperature !temperature
  DOUBLEPRECISION, DIMENSION(N**2) :: distance
  
  seed = 918172
  CALL SRAND(seed)

  OPEN(unit=11,file="position.dat")
  OPEN(unit=12,file="force.dat")
  OPEN(unit=14,file="accel.dat")
  OPEN(unit=15,file="energy.dat")
  OPEN(unit=16,file="vcm.dat")
  OPEN(unit=17,file="temperature.dat")
  OPEN(unit=18,file="velocity.dat")
  OPEN(unit=19,file="distance.dat")
  OPEN(unit=20,file="vx.dat")
!=============set initial position and velocities====================================
  DO i = 1, N
   !PRINT*, 0.05D0*(2D0*RAND()-1D0)
    px(i) =  0.5D0 + MOD(i-1,10) * sidelen/DBLE(10) + 0.05D0*(2D0*RAND()-1D0)
    py(i) =  0.5D0 + INT((i-0.1)/10D0) * sidelen/DBLE(10) + 0.05D0*(2D0*RAND()-1D0)
	vx(i) = v0
	vy(i) = v0
!    WRITE(11,*) i, px(i), py(i)
  ENDDO

!Euler to get px_old
 DO i = 1, N
   px_older(i) = px(i)
   py_older(i) = py(i)
 ENDDO
 DO i = 1, N
   px_old(i) = px_older(i) + vx(i)*deltaT
   py_old(i) = py_older(i) + vy(i)*deltaT
 ENDDO
 

!Verlet
  DO it = 1, maxT
    CALL forceij

    DO i = 1, N
	  CALL accel(i,ax(i),ay(i))
	  WRITE(14,*) i, ax(i), ay(i)
	  px(i) = 2D0*px_old(i) - px_older(i) + ax(i)*deltaT**2
	  IF ( DABS(px(i)) .GT. 2D0*sidelen ) PRINT*, "TROUBLE: FOUL BALL!"
!	  IF ( px(i) .GT. sidelen ) THEN
!	    px(i) = px(i) - sidelen
!	  ELSEIF ( px(i) .LT. 0D0 ) THEN
!	    px(i) = px(i) + sidelen
!	  ELSE                              
!	    px(i) = px(i)
!      ENDIF

	  py(i) = 2D0*py_old(i) - py_older(i) + ay(i)*deltaT**2
	  IF ( DABS(py(i)) .GT. 2D0*sidelen ) PRINT*, "TROUBLE: FOUL BALL!"
!	  IF ( py(i) .GT. sidelen ) THEN
!	    py(i) = py(i) - sidelen
!	  ELSEIF ( py(i) .LT. 0D0 ) THEN
!	    py(i) = py(i) + sidelen
!	  ELSE                              
!	    py(i) = py(i)
!      ENDIF
	  vx(i) = ( px(i) - px_older(i) )/( 2D0 * deltaT )
	  vy(i) = ( py(i) - py_older(i) )/( 2D0 * deltaT )
IF (vx(i) > 100) PRINT*, it*deltaT, i, vx(i)
      px_older(i) = px_old(i)
	  py_older(i) = py_old(i)
	  px_old(i) = px(i)
	  py_old(i) = py(i)

    ENDDO !i

    !calculating energy, vcm, temperature
	ke = 0D0
	vcmx = 0D0
	vcmy = 0D0
	DO i = 1, N
      ke = ke + 0.5D0 * ( (vx(i))**2 + (vy(i))**2 )
	  vcmx = vcmx + vx(i)
	  vcmy = vcmy + vy(i)
	ENDDO
	temperature = ke/DBLE(N) !actually kT
!	WRITE(*, *) "KE = ", ke, ", PE = ", pe, ", PE+KE = ", pe+ke
	WRITE(15,'(4ES20.10)') it*deltaT, ke, pe, pe+ke
	WRITE(16,'(4ES20.10)') it*deltaT, vcmx, vcmy, DSQRT(vcmx**2+vcmy**2)
	WRITE(17,'(2ES20.10)') it*deltaT, temperature 

  ENDDO !it

!record position
  DO i = 1, N
    WRITE(11,*) i, px(i), py(i)
	WRITE(20,*) i, vx(i)
  ENDDO

!record distance
  CALL rcorrelation

  PRINT*, "Message: Position of each particle recorded in file ""position.dat"", in order of ""label, x, y"" "
  PRINT*, "Message: Force of each pair of particles recorded in file ""force.dat"", in order of ""i, j, forcex, forcey"" "
  PRINT*, "Message: Acceleration of each particles recorded in fiel ""accel.dat"", in order of ""label, ax, ay"" "
  PRINT*, "Message: Energies of the system recorded in ""energy.dat"", in order of ""time, kinetic, potential and total energy"". "
  PRINT*, "Message: Center of mass velocities recorded in ""vcm.dat"", in order of ""time, x comp, y comp, and magnitude of CM velocity."" "
  PRINT*, "Message: Temperature recorded in ""temperature.dat"", in order of ""time, temperature."" "
  PRINT*, "Message: Final vx recorded in ""vx.dat"", in order of ""label, vx."" "


! PAUSE

 CONTAINS

 SUBROUTINE forceij     !calculate Fij, the force of i ON j, of each configuration and store them into the array
   IMPLICIT NONE
   INTEGER :: ifor, jfor

   DOUBLEPRECISION, PARAMETER :: rcut = 3D0

   DOUBLEPRECISION :: x, y  !For calculating force
   DOUBLEPRECISION :: r
   DOUBLEPRECISION :: costhe, sinthe  !cos and sin of theta

   pe = 0D0
   DO ifor = 1, N
     DO jfor = 1, ifor      !so jfor <= ifor
       IF ( ifor == jfor ) THEN
         forcex(ifor,jfor) = 0D0
		 forcey(ifor,jfor) = 0D0
       ELSE 
         IF ( (px_old(jfor) - px_old(ifor)) > (sidelen/2D0) ) THEN       !right -> left
           x = px_old(jfor) - sidelen
         ELSEIF ( (px_old(jfor) - px_old(ifor)) < -sidelen/2D0 ) THEN  !left -> right
           x = px_old(jfor) + sidelen
         ELSE                                !no telegraph
           x = px_old(jfor)
         ENDIF

         IF ( (py_old(jfor) - py_old(ifor)) > sidelen/2D0 ) THEN       !above -> under
           y = py_old(jfor) - sidelen
         ELSEIF ( (py_old(jfor) - py_old(ifor)) < -sidelen/2D0 ) THEN  !under -> above
           y = py_old(jfor) + sidelen
         ELSE                                !no telegraph
           y = py_old(jfor)
         ENDIF

         r = distBetween(px_old(ifor),py_old(ifor),x,y)
         IF ( r > rcut ) THEN
           forcex(ifor,jfor) = 0D0
		   forcey(ifor,jfor) = 0D0
         ELSE
           costhe = (x - px_old(ifor))/r   !the force by j on particle i
           sinthe = (y - py_old(ifor))/r
           forcex(ifor,jfor) = ljForce(r)*costhe
		   forcey(ifor,jfor) = ljForce(r)*sinthe
		 ENDIF !rcut
	   ENDIF !(ifor ==jfor)
	   pe = pe + ljPotential(r)
     ENDDO !jfor
   ENDDO !ifor  

   DO ifor = 1, N
     DO jfor = ifor, N      !so jfor >= ifor
       forcex(ifor,jfor) = -forcex(jfor,ifor)
	   forcey(ifor,jfor) = -forcey(jfor,ifor)
     ENDDO !jfor
   ENDDO !ifor 
 END SUBROUTINE

 SUBROUTINE accel(i,ax,ay)  !inquiring the i_th particle, ax and ay are the outputs
   IMPLICIT NONE
   DOUBLEPRECISION, PARAMETER :: m = 1D0
   DOUBLEPRECISION :: ax, ay
   INTEGER :: i,j
   ax = 0D0
   ay = 0D0
   DO j = 1, N
!     PRINT*, forcex(j,i)
     ax = ax + forcex(j,i)/m
	 ay = ay + forcey(j,i)/m
   ENDDO
!   PRINT*, ax, ay
   RETURN
 END SUBROUTINE

   
 DOUBLEPRECISION FUNCTION distBetween(x1,y1,x2,y2)
   IMPLICIT NONE
   DOUBLEPRECISION :: x1,y1,x2,y2
   distBetween = DSQRT( (x1-x2)**2 + (y1-y2)**2 )
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION ljForce(r)
   IMPLICIT NONE
   DOUBLEPRECISION :: r
   ljForce = 24D0 * ( 2D0 / r**13 - 1D0 / r**7)
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION ljpotential(r)
   IMPLICIT NONE
   DOUBLEPRECISION :: r
   ljpotential = 4D0*( 1D0 / r**12 - 1D0 / r**6 )
   RETURN
 END FUNCTION


 SUBROUTINE rcorrelation     !recording rij
   IMPLICIT NONE
   INTEGER :: ifor, jfor

   DOUBLEPRECISION :: x, y  !For calculating force
   DOUBLEPRECISION :: r

   DO ifor = 1, N
     DO jfor = 1, ifor      !so jfor <= ifor
       IF ( ifor == jfor ) THEN
	     CONTINUE
       ELSE 
         IF ( (px_old(jfor) - px_old(ifor)) > (sidelen/2D0) ) THEN       !right -> left
           x = px_old(jfor) - sidelen
         ELSEIF ( (px_old(jfor) - px_old(ifor)) < -sidelen/2D0 ) THEN  !left -> right
           x = px_old(jfor) + sidelen
         ELSE                                !no telegraph
           x = px_old(jfor)
         ENDIF

         IF ( (py_old(jfor) - py_old(ifor)) > sidelen/2D0 ) THEN       !above -> under
           y = py_old(jfor) - sidelen
         ELSEIF ( (py_old(jfor) - py_old(ifor)) < -sidelen/2D0 ) THEN  !under -> above
           y = py_old(jfor) + sidelen
         ELSE                                !no telegraph
           y = py_old(jfor)
         ENDIF

         r = distBetween(px_old(ifor),py_old(ifor),x,y)
		 WRITE(19,*) r
	   ENDIF !(ifor ==jfor)
     ENDDO !jfor
   ENDDO !ifor  
 END SUBROUTINE


END PROGRAM



