program main
 IMPLICIT NONE
 DOUBLEPRECISION :: richard
 INTEGER :: k, nlow, ncap


 ncap = 2!ncap is the number of Pade temrs being grouped together to form Richardson transformation. The N index in the expression of Q0 
 WRITE(1,*) "N for Richardson is", ncap
 WRITE(1,*) "n,       Pade,       Richardson,       exp(-10)"

do nlow=1,20 !nlow is the N index of P_NN
!do Richardson series
   richard = 0D0
  DO k = 0, ncap 
   richard = richard + pade(2*(k+nlow)+1) *dble(nlow+k)**ncap * (-1D0)**(k+ncap)/dble( factorial(k)*factorial(ncap-k) )
  ENDDO 

  WRITE(1, '(I2, 1F20.15, 2ES25.14)') nlow, pade(2*nlow+1), richard, DEXP(-10D0)
end do


contains

double precision function pade(pn)
integer,parameter:: NUM=60
integer::n,m,i,flag,pn
double precision::x=10d0,Avalue,Bvalue,Avalue_1,Avalue_2,Bvalue_1,Bvalue_2
double precision, dimension(0:NUM-1):: e,q
double precision, dimension(1:NUM):: a

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
end do
a(flag)=-q(0)

!calculate es and qs and get the value of as
do i=3,NUM,2
   flag=i
   do n=0,NUM-flag
      e(n)=q(n+1)-q(n)+e(n+1) 
   end do
   a(flag)=-e(0)
  
   m=m+1
   flag=flag+1
   do n=0,NUM-flag
      q(n)=e(n+1)/e(n)*q(n+1)
   end do  
   a(flag)=-q(0)
end do

!STEP 2, evaluate using pade
   Avalue_1=0d0
   Avalue_2=1d0
   Bvalue_1=1d0
   Bvalue_2=0d0  
   do i=1,pn
      if (i==1) then
         Avalue=Avalue_1+a(i)*Avalue_2
         Bvalue=Bvalue_1+a(i)*Bvalue_2
      else
         Avalue=Avalue_1+a(i)*x*Avalue_2
         Bvalue=Bvalue_1+a(i)*x*Bvalue_2
      end if
      Avalue_2=Avalue_1
      Avalue_1=Avalue
      Bvalue_2=Bvalue_1
      Bvalue_1=Bvalue
   end do
   pade=Avalue/Bvalue
end function pade


integer FUNCTION factorial(nn)
  IMPLICIT NONE
  INTEGER :: iter,nn,fact

  fact=1
  DO iter = 1, nn
   fact = fact*DBLE(iter)
  ENDDO
   factorial = fact
   RETURN
 END FUNCTION

 
end program
