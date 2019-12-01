program limits
  IMPLICIT NONE
  real:: realUnder, realOver !parameters for single precision limit
  double precision:: doubleUnder,doubleOver !parameters for double precision limit 
 integer:: N,loopNum
  OPEN(unit=2,file='data2.txt')
  OPEN(unit=3,file='data3.txt') 
 realUnder=1.0
  realOver=1.0
  doubleUnder=1.0
  doubleOver=1.0
  N=2000
  do loopNum=1, N
     realUnder=realUnder/2.0 
     realOver=realOver*2.0    
     doubleUnder=doubleUnder/2.
     doubleOver=doubleOver*2.
     write (2,*), loopNum, realUnder, realOver
     write (3,*), loopNum, doubleUnder, doubleOver
  end do
  close(2)
  close(3)
end program limits
