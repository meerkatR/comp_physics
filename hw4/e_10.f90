!     Purpose:
!     Using Pade P_{NN} to approximate exp(-10) with array construct and compare with Richardson
!  
!
!
!     Record of revisions:
!     Date             Programmer           Description of change
!     ===========      ==========           =====================
!     5.Oct.2009      Jen Hsu            

PROGRAM padeapprox
 IMPLICIT NONE
 INTEGER :: niter
 DOUBLEPRECISION :: xmax = 50D0
 DOUBLEPRECISION :: x, dx = 0.1D0
 DOUBLEPRECISION :: relerror
 DOUBLEPRECISION :: richsum
 INTEGER :: k, nlow, ncap

!vary N
 nlow = 5
 OPEN(unit=19, file='exp_10_varyNcap.dat')
 WRITE(19,*) "n is fixed and N is varied."
 WRITE(19,*) "n, the index of the beginning of Richardson transform, Pnn, is", nlow
 WRITE(19,*) "N,       Pade,       Richardson,       exp(-10),          relative error"

 DO niter = 1, 30  !this is index in P_NN
  relerror = ( expAppro(2*niter+1,10D0) - DEXP(-10D0) )/DEXP(-10D0)
  
!do Richardson series
  richsum = 0D0
  DO k = 0, niter 
   richsum = richsum + expAppro(2*(k+nlow)+1, 10D0) * (nlow+k)**niter * (-1D0)**(k+niter)/( factorial(k)*factorial(niter-k) )
!   aSum = aSum + harsum(nlow+k) * (nlow+k)**niter * (-1D0)**(k+niter)/( factorial(k)*factorial(niter-k) )
  ENDDO 

  WRITE(19, '(I2, 1F20.15, 3ES25.14)') niter, expAppro(2*niter+1,10D0), richsum, DEXP(-10D0), DABS(relerror)
 ENDDO


!varying n
 ncap = 2
 OPEN(unit=18, file='exp_10_varynlow.dat')
 WRITE(18,*) "N is fixed and n is varied."
 WRITE(18,*) "N, the number of Pade terms being grouped together to form Richardson transform, is", ncap
 WRITE(18,*) "n,       Pade,       Richardson,       exp(-10),          relative error"

 DO nlow = 1, 20  !this is index in P_NN
  relerror = ( expAppro(2*nlow+1,10D0) - DEXP(-10D0) )/DEXP(-10D0)

!do Richardson series
  richsum = 0D0
  DO k = 0, ncap 
   richsum = richsum + expAppro(2*(k+nlow)+1, 10D0) * (nlow+k)**ncap * (-1D0)**(k+ncap)/( factorial(k)*factorial(ncap-k) )
  ENDDO 

  WRITE(18, '(I2, 1F20.15, 3ES25.14)') nlow, expAppro(2*nlow+1,10D0), richsum, DEXP(-10D0), DABS(relerror)
 ENDDO



 PRINT*, "Message: Approx result written in 'exp_10_varyN.dat' file, in format of 'N', 'expAppro(P_NN)', 'Richardson to nlow',  &
& 'exp(-10)(ans)', 'error' "
 PRINT*, "Message: Approx result written in 'exp_10_varyn.dat' file, in format of 'n', 'expAppro(P_nn)', 'Richardson to nlow',  &
& 'exp(-10)(ans)', 'error' "



 PAUSE

 CONTAINS
 DOUBLEPRECISION FUNCTION expAppro(N,x)
 INTEGER :: niter, miter, aiter
 INTEGER :: N      !a up to a(10) corresponds to P_{N,N}
 INTEGER :: Nquo

 DOUBLEPRECISION :: x
 DOUBLEPRECISION, PARAMETER :: pi = 3.14159265358979323846264338328
 DOUBLEPRECISION :: result, relError
 DOUBLEPRECISION, DIMENSION(0:N+1,0:N+1) :: e, q
 DOUBLEPRECISION, DIMENSION(0:N+1) :: taylorc, padea !co-eff of Taylor, co-eff of Pade
 DOUBLEPRECISION, DIMENSION(-1:N) :: acap, bcap
 DOUBLEPRECISION, DIMENSION(N) :: a     

 Nquo = N

 DO miter = 0, Nquo+1
  DO niter = 0, Nquo+1
   e(miter,niter) = 0D0
   q(miter,niter) = 0D0
  ENDDO
 ENDDO

 DO miter = 0, Nquo
  DO niter = Nquo, 0, -1
   IF (miter == 0) THEN
    e(miter,niter) = 0D0
   ELSE IF (miter == 1) THEN
    q(miter,niter) = -1D0/DBLE(niter+1)
    e(miter,niter) = q(miter,niter+1) - q(miter,niter) + e(miter-1,niter+1)
   ELSE
	q(miter,niter) = e(miter-1, niter+1)/e(miter-1,niter)*q(miter-1,niter+1)    
    e(miter,niter) = q(miter,niter+1) - q(miter,niter) + e(miter-1,niter+1)
!	q(miter+1,niter) = e(miter, niter+1)/e(miter,niter)*q(miter,niter+1)
   ENDIF
  ENDDO !niter
 ENDDO !miter

 a(1) = 1D0
 DO aiter = 3, N, 2
  a(aiter) = -e(NINT(DBLE(aiter-1)/2D0),0)*x
 ENDDO
 DO aiter = 2, N, 2
  a(aiter) = -q(NINT(DBLE(aiter)/2D0),0)*x
 ENDDO

  !transform a and b's into A and B's
  acap(-1) = 1D0
  acap(0)  = 0D0
  bcap(-1) = 0D0
  bcap(0) =  1D0
  DO aiter = 1, N
   acap(aiter) = acap(aiter-1)+a(aiter)*acap(aiter-2)
   bcap(aiter) = bcap(aiter-1)+a(aiter)*bcap(aiter-2)
  ENDDO

  expAppro = acap(N)/bcap(N)
  RETURN
  END FUNCTION

 DOUBLEPRECISION FUNCTION factorial(nn)
  IMPLICIT NONE
  INTEGER :: iter,nn
  DOUBLEPRECISION :: fact

  fact=1D0
  DO iter = 1, nn
   fact = fact*DBLE(iter)
  ENDDO
   factorial = fact
   RETURN
 END FUNCTION

 DOUBLEPRECISION FUNCTION harsum(nn)
  IMPLICIT NONE
  INTEGER :: iter, nn
  DOUBLEPRECISION :: sum

  sum = 0D0
  DO iter = 1, nn
   sum = sum + 1D0/iter**2
  ENDDO
  harsum = sum
  RETURN
 END FUNCTION
 
 DOUBLEPRECISION FUNCTION expsum(nn)
  IMPLICIT NONE
  INTEGER :: iter, nn
  DOUBLEPRECISION :: sum

  sum = 1D0
  DO iter = 1, nn
   sum = sum + (-10D0)**iter/factorial(iter)
  ENDDO
  expsum = sum
  RETURN
 END FUNCTION
   
END PROGRAM
