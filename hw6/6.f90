program ocean
IMPLICIT NONE
integer,parameter:: N=1000,M=320
integer::i,j
double precision::ul=15d0, ur=20d0
double precision:: delx, delt,delv
double precision::alpha=8d-7, beta=207d-6
double precision,dimension(-N:N):: u,uold

delx=1d0/N
delt=0.5d0
!delt=0.630 !test of stability 
!delx=0.0009d0!test of stability

!initial T for sea
do i=-N,0
   uold(i)=ul
   u(i)=ul
end do
!initial T for sky
do i=1,N
   uold(i)=ur
   u(i)=ur
end do

delv=0d0
do j=1,M !loop of time 
   do i=-(N-1),N-1 !loop of space
      u(i)=uold(i)+delt*alpha*((uold(i+1)-2*uold(i)+uold(i-1))/delx**2)
      !calculate deltaV for each disc for sea
      if (i<=0) then
         delv=delv+delx*beta*(u(i)-uold(i))
      end if
   end do
   write(2,*) j*delt,delv

   !store old value of Temprature distribution
   do i=-N,N
      uold(i)=u(i)
   end do
end do

delv=0d0
write(1,*) "Temprature distribution @ t=",M*delt
do i=-N,N
   write(1,*) i*delx,u(i) 
end do

end program

