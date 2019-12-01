program integral
IMPLICIT NONE
double precision:: a,b,h,analvalue,trapezoid,simpson,gauss,trapErr,simpErr,gausErr 
double precision, allocatable, dimension(:) :: w,x
double precision:: fcn
integer:: N,i


a=0d0
b=1.d0
analValue=1.d0-exp(-1.d0)

do N=5, 400, 4
   ALLOCATE(w(N))
   ALLOCATE(x(N))
   h=(b-a)/(N-1)

!Trapezoid Method
  w(1)=h/2.d0
  w(N)=h/2.d0
  do i=2, N-1
     w(i)=h
  end do

  do i=1, N
     x(i)=a+(i-1)*h
  end do

  trapezoid=0.d0
  do i=1,N
     trapezoid=trapezoid+w(i)*fcn(x(i))
  end do
  trapErr=DABS(trapezoid-analValue)/analValue


!Simpson's rule
  w(1)=h/3.d0
  w(N)=h/3.d0
  do i=2, N-1, 2
     w(i)=4.d0*h/3.d0
  end do
 
  do i=3, N-1, 2
     w(i)=2.d0*h/3.d0
  end do

  simpson=0.d0
  do i=1,N
     simpson=simpson+w(i)*fcn(x(i))
  end do
  simpErr=DABS(simpson-analValue)/analValue

!Gaussioan quadrature methods
  CALL GaussLegendre(x,w,N,a,b)
  gauss=0.d0
  do i=1,N
     gauss=gauss+w(i)*fcn(x(i))
  end do
  gausErr=DABS(gauss-analValue)/analValue
       
  IF(ALLOCATED(w)) DEALLOCATE(w)
  IF(ALLOCATED(x)) DEALLOCATE(x)
  write(1,*) N,trapErr,simpErr,gausErr
end do
end program 
INCLUDE 'GaussLegendre.f'


double precision function fcn(x)
IMPLICIT NONE
double precision:: x
fcn=exp(-x)
return
end function fcn
