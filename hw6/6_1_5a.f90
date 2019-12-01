program SORJacobi
IMPLICIT NONE
integer, parameter::N=40
integer:: i,j,l,neg,pos,NUM
double precision:: delta,norm,vstar,sum2
double precision:: PI,rjac,w
double precision, dimension(-N:N, -N:N):: v
double precision, dimension(-N:N):: prev, curr !in order to perform Jacobi Method, record the previous row valude and current row value

PI = 4.d0*datan(1.d0)

delta=1d0/N  !grid spacing
neg=nint(-0.25d0/delta) !the grid position of left plate
pos=nint(0.25d0/delta)  !the grid position of right plate

!initial distribution of v
do j=-N, N
   do l=-N, N
      if (j==neg .AND. l>=neg .AND. l<=pos) then
         v(j,l)=-1d0
      else if (j==pos .AND. l>=neg .AND. l<=pos) then 
         v(j,l)=1d0
      else 
         v(j,l)=0d0 
      end if
   end do
end do


do i=-N, N
   prev(i)=0d0 !the initial value of prev array=v(j,-N)=0
   curr(i)=0d0 !main purpose is to set the boundary value for curr array
end do

norm=1d0
NUM=0

rjac=1d0-PI**2/(2d0*dble(2*N+1)**2)
w=2d0/(1+dsqrt(1-rjac**2))

do while(norm>=1d-8 .AND. NUM<10000) !Add a new looping condition in case the method does't converge
   NUM=NUM+1 !iteration number counting
   sum2=0d0
   do l=-(N-1), (N-1) !sweep along rows
      do j=-(N-1), (N-1)
         curr(j)=v(j,l) !use curr array to record old v(j,l)
         if (j==neg .AND. l>=neg .AND. l<=pos) then
            continue
         else if (j==pos .AND. l>=neg .AND. l<=pos) then 
            continue
         else
           vstar=0.25d0*(v(j+1,l)+curr(j-1)+v(j,l+1)+prev(j))
           v(j,l)=w*vstar+(1d0-w)*curr(j)
         end if
         sum2=sum2+(v(j,l)-curr(j))**2!square sum of difference matrix
      end do !j loop
      do i=-N,N
         prev(i)=curr(i) !use prev array to record the old row
      end do
    end do !l loop
    norm=dsqrt(sum2)!Frobenius norm of difference matrix
    write(1,*) NUM,norm
end do !do while loop
print *,NUM


end program 
