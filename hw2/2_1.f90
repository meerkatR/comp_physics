program precision
real:: eps, one !for double precision, change 'real' to 'real(8)'
integer:: i, N
eps=1.
N=100
OPEN(unit=1,file='data1.txt')
do i=1, N
   eps=eps/2.
   one=1.+eps
   write(1,'(1X,I5,F28.23,ES30.6)') i,one,eps
end do
close(1)
end program precision
