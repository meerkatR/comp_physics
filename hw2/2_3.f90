program main
IMPLICIT NONE
integer:: n,dFact
integer:: i,M
double precision::x,dx,jRecurs,jAsympt,sphBessel
double precision::error,relativeError

open(unit=1,file='data4_1.txt')
open(unit=2,file='data4_2.txt')
print *,'Enter n for Spherical Bessel functions, jn(x):'
read *, n

x=0.0
M=500
dx=(1.0-0.0)/M

do i=1,M
   x=x+dx
   jRecurs=sphBessel(n,x)
   jAsympt=x**n/dFact(2*n+1)*(1-x**2/(2*(2*n+3)))
   error=jRecurs-jAsympt
   relativeError=error/jAsympt
   write(1,*), x,jRecurs,jAsympt
   write(2,*), x,relativeError
end do

close(1)
close(2)

end program main




recursive integer function dFact(n) result(Fvalue)
  IMPLICIT NONE
  integer:: n, Fvalue

  if (n==0) then
    Fvalue=1
  else if (n==1) then
    Fvalue=1
  else
    Fvalue=dFact(n-2)*n
  end if

  return
end function dFact




recursive double precision function sphBessel(n,x) result(Bvalue)
  IMPLICIT NONE
  integer:: n
  double precision:: x,Bvalue

  if (n.eq.0) then
    Bvalue=dsin(x)/x
  else if (n.eq.1) then
    Bvalue=dsin(x)/(x**2)-dcos(x)/x
  else
    Bvalue=(2*(n-1)+1)/x*sphBessel(n-1,x)-sphBessel(n-2,x)
  end if
  return
end function sphBessel
