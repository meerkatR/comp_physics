      subroutine GaussLegendre(x,w,n,lowlim,uplim)
!c     **********************************************************
!c     *                                                        *
!c     *   Calculate weights and abscissae for Gauss-Legendre   *
!c     *   quadrature.                                          *
!c     *   These are to be used in the formula                  *
!c     *                                                        *
!c     *     Intgl_{lowlim}^{uplim}  f(x) dx                    *
!c     *          =   Sum_{j=1 to N} w(j)*f(x_j)                *
!c     *                                                        *
!c     **********************************************************

      Implicit none
      integer n
      real*8 lowlim,uplim, w(n),x(n)
      real*8 eps
      parameter(eps=3.0d-14)

      integer i,j,m
      real*8 p1,p2,p3,pp,xl,xm,z,z1

      m= (n+1)/2
      xm= 0.5d0*(uplim + lowlim)
      xl= 0.5d0*(uplim - lowlim)
      do i=1,m
        z= dcos(3.141592653589793d0*(i-0.25d0)/(0.5d0+n))
1       continue 
        p1= 1.0d0
        p2= 0.0d0
        do j=1,n
          p3= p2
          p2= p1
          p1= ( (2.0d0*j -1.0d0)*z*p2 +(1.0d0 -j)*p3 )/j
          end do
        pp= n*(z*p1 -p2)/(z*z -1.0d0)
        z1= z
        z= z1 -p1/pp

        if (dabs(z-z1) .gt. eps) goto 1
        x(i)= xm -xl*z
        x(n+1-i)= xm +xl*z
        w(i)= 2.0d0*xl/( (1.0d0-z*z)*pp*pp )
        w(n+1-i)= w(i)
        end do
      return
      end
