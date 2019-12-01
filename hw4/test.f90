program main
double precision::a,b
print *, a(b,5.0)
end program main





double precision function a(b,x) 
IMPLICIT NONE
double precision::x
double precision, EXTERNAL::b

a=b(x+1)/b(x)
end function

double precision function b(x)
IMPLICIT NONE
double precision::x
!double precision, EXTERNAL::a

b=x
end function

