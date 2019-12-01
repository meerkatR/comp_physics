program
IMPLICIT NONE
integer, parameter::seed=918172
integer,parameter::N=100,step=1000,R=10
double precision::RAND
double precision::v0=0d0,dt=0.005d0,deltar=0.05d0
double precision::forcex,forcey
double precision,dimension(N)::gx(i),gy(i),xnew,xcurr,xprev,ynew,ycurr,yprev
double precision,dimension(N)::vx,vy,ax,ay

!initializing
do i=1,R
   gx(i)=0.5+(R-1)
   gy(i)=0.5+(R-1)
end do
call srand(seed)
do i=1,N
   xcurr(i)=gx(i)+2d0*(rand()-0.5d0)*deltar
   ycurr(i)=gy(i)+2d0*(rand()-0.5d0)*deltar
   vx(i)=v0
   vy(i)=v0
   !we need xprev, yprev to run Verlet
   xprev(i)=xcurr(i)-vx(i)*dt
   yprev(i)=ycurr(i)-vy(i)*dt
end do

do i=1,N
   force=0d0
   do j=1,N
      if (j.LE.i) then
         !calculate r
         !calculate theta

         if(r>rcut) then
           f(i,j)=0
         else
          ! force
           f(i,j)=24d0*(2d0/r**13-1d0/r**7)
          
          ax=force*cos_theta
          ay=force*sin_theta
         end if
      end if   
   end do!all other particles

  xnew(i)=2d0*xcurr(i)-xprev(i)+ax
  ynew(i)=2d0*ycurr(i)-yprev(i)+ay

      
end do!go through each particle




end program
