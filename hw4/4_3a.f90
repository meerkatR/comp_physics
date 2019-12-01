program quotientDifference
integer,parameter:: NUM=10
integer::n,m,i,flag,pn,step
double precision::x,xmax,dx,pade,pA,pB
double precision, dimension(0:NUM-1):: e,q
double precision, dimension(1:NUM):: a
common /data1/ a

a(1)=1d0

flag=1
m=0
e=0d0

flag=flag+1
m=m+1
do n=0,NUM-flag
   q(n)=-1d0/(dble(n)+1d0)
print *, m,n, q(n)
end do
a(flag)=-q(0)


do i=3,NUM,2
   flag=i
   do n=0,NUM-flag
      e(n)=q(n+1)-q(n)+e(n+1) 
  ! print *, M,n, 'e=',e(n)
   end do
   a(flag)=-e(0)
  
   m=m+1
   flag=flag+1
   do n=0,NUM-flag
      q(n)=e(n+1)/e(n)*q(n+1)
   !print *, m,n, 'q=',q(n)
!   print *, flag, a(flag)
   end do  
   a(flag)=-q(0)
end do

do i=1, NUM
   print*, a(i)
end do

print *, 'enter the value of N for Pade:'
read *, pn

x=0d0
xmax=50d0
step=1000
dx=(xmax-0d0)/step

do i=0,step
   pade=pA(pn,x)/pB(pn,x)
   write(1,*) x,dexp(-x),pade
   x=x+dx
end do

end program quotientDifference


recursive double precision function pA(pn,x) result(Avalue)
IMPLICIT NONE
integer, parameter:: NUM=10
integer:: pn
double precision:: x, Avalue
double precision, dimension(1:NUM)::a
COMMON /data1/ a

if (pn==-1) then
   Avalue=1d0
else if (pn==0) then
   Avalue=0D0
else if (pn==1) then
   Avalue=pA(pn-1,x)+a(pn)*pA(pn-2,x)
else
   Avalue=pA(pn-1,x)+a(pn)*x*pA(pn-2,x)
end if
return
end function pA

recursive double precision function pB(pn,x) result(Bvalue)
IMPLICIT NONE
integer, parameter::NUM=10
integer:: pn
double precision:: x, Bvalue
double precision,dimension(1:NUM)::a
COMMON /data1/ a

if (pn==-1) then
   Bvalue=0D0
else if (pn==0) then
   Bvalue=1D0
else if(pn==1) then
   Bvalue=pB(pn-1,x)+a(pn)*pB(pn-2,x)
else
   Bvalue=pB(pn-1,x)+a(pn)*x*pB(pn-2,x)
end if
return
end function pB
 
