program quotientDifference
integer:: m,n
double precision::res, e,q
!do m=0, 10
!   do  n=0, 10

call e(q,0,0)
res=e
print *, res
!   print *, e(q,0,0)
!   end do
!end do
end program quotientDifference 



recursive double precision function e(q,m,n) result(evalue)
IMPLICIT NONE
double precision, EXTERNAL:: q
integer:: m, n
double precision:: evalue
if (m==0) then 
   evalue=0d0
else
   evalue=q(e,m,n+1)-q(e,m,n)+e(q,m-1,n-1)
end if
return
end function e

recursive double precision function q(e,m,n) result(qvalue)
IMPLICIT NONE
double precision, EXTERNAL:: e
integer:: m,n
double precision:: qvalue
double precision, dimension(0:9):: c

if (m==1) then
   qvalue=c(n+1)/c(n)
else
   qvalue=e(q,m,n+1)/e(q,m,n)*q(e,m,n+1)
end if
return
end function q


