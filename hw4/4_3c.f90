program quotientDifference
integer,parameter:: NUM=20
integer::n,m,i,j,flag,pn,step
double precision::x,xmax,dx,error
double precision::pade,Avalue,Bvalue,Avalue_1,Avalue_2,Bvalue_1,Bvalue_2
double precision, dimension(0:NUM):: e,q
double precision, dimension(1:NUM+1):: a


!STEP 1, get an for cotinued fraction using quotient-difference algorithm
a(1)=1d0

!initial value for es
flag=1!flag is the index of array a, and also the total column of e and q
m=0!m is the lower index of e and p
e=0d0

!initial value for qs
flag=flag+1
m=m+1
do n=0,NUM-flag
   q(n)=-1d0/(dble(n)+1d0)
print *, m,n,'q=', q(n)
end do
a(flag)=-q(0)

!calculate es and qs and get the value of as
do i=3,NUM,2
   flag=i
   do n=0,NUM-flag
      e(n)=q(n+1)-q(n)+e(n+1) 
   print *, m,n, 'e=',e(n)
   end do
   a(flag)=-e(0)
  
   m=m+1
   flag=flag+1
   do n=0,NUM-flag
      q(n)=e(n+1)/e(n)*q(n+1)
   print *, m,n, 'q=',q(n)
   end do  
   a(flag)=-q(0)
end do

do i=1, NUM
   print*, a(i)
end do


!STEP 2, evaluate using pade
print *, 'enter the value of N for Pade:'
read *, pn

x=0d0
xmax=40d0
step=1000
dx=(xmax-0d0)/step

do i=0,step
   Avalue_1=0d0
   Avalue_2=1d0
   Bvalue_1=1d0
   Bvalue_2=0d0  
   do j=1,pn
      if (j==1) then
         Avalue=Avalue_1+a(j)*Avalue_2
         Bvalue=Bvalue_1+a(j)*Bvalue_2
      else
         Avalue=Avalue_1+a(j)*x*Avalue_2
         Bvalue=Bvalue_1+a(j)*x*Bvalue_2
      end if
      Avalue_2=Avalue_1
      Avalue_1=Avalue
      Bvalue_2=Bvalue_1
      Bvalue_1=Bvalue
   end do
   pade=Avalue/Bvalue
   error=(pade-dexp(-x))/dexp(-x)
   write(1,*) x,dexp(-x),pade,error
   x=x+dx
end do

end program quotientDifference
