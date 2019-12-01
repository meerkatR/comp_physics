program sphereVolume
IMPLICIT NONE
integer, parameter:: M=10000
integer:: seed, N, numIn
integer:: i, numTimes, numThrow
real(8):: x, y, z, r 
real(8):: r0, cube, sphere, rm, sRange 
real(8), dimension(M):: volume

open(unit=1, file='data1.txt') !if N=100,000, file='data2.txt'

volume=0.0

N=1000
r0=1.0


sRange=2*r0
rm=r0-sRange
cube=sRange**3

seed=918172
call srand(seed)


do numTimes=1, M
   numIn=0  
   do numThrow=1, N
      x=(rand()*sRange)+rm
      y=(rand()*sRange)+rm
      z=(rand()*sRange)+rm
      r=sqrt(x**2+y**2+z**2)
      if (r<=r0) then
         numIn=numIn+1
      end if
   end do
   sphere=real(numIn)/N*cube
   volume(numTimes)=sphere
 write(1,*) volume(numTimes)
end do

close(1)
end program sphereVolume




