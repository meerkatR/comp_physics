program GaussSeidel
IMPLICIT NONE
integer, parameter::N=80
integer:: i,j,k,l,neg,pos,NUM
double precision:: delta,old,sum2,norm
double precision, dimension(-N:N, -N:N):: v

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

norm=1d0
NUM=0

do while(norm>=1d-8)
   NUM=NUM+1 !iteration number counting
   sum2=0d0
   do l=-(N-1), (N-1) !sweep along rows
      do j=-(N-1), (N-1)
         old=v(j,l) !record the old value of v(j,l)
         if (j==neg .AND. l>=neg .AND. l<=pos) then
            continue
         else if (j==pos .AND. l>=neg .AND. l<=pos) then 
            continue
         else
           v(j,l)=0.25d0*(v(j+1,l)+v(j-1,l)+v(j,l+1)+v(j,l-1))
         end if
         sum2=sum2+(v(j,l)-old)**2!square sum of difference matrix
      end do !j loop
    end do !l loop
    norm=dsqrt(sum2)!Frobenius norm of difference matrix
    write(1,*) NUM,norm
end do !do while loop
print *,NUM

end program 
