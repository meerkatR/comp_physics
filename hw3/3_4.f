c     vegasEx.f    
      program vegasIntegral
      implicit none
      real*8 reg(12),pi,result,analyresult,sdr,chi
      real*8 negInf, posInf,analysresult,error
      integer ndim,ncall,itmx,nprn,idum
      common /ranno/ idum
      external fcn

      pi = 4.d0*datan(1.d0)

      idum=9812371
      ndim=6
      itmx=10
      nprn=-1

      negInf=-2.0
      posInf=2.0

      reg(1) = negInf
      reg(2) = negInf
      reg(3) = negInf
      reg(4) = negInf
      reg(5) = negInf
      reg(6) = negInf
      reg(7) = posInf
      reg(8) = posInf
      reg(9) = posInf
      reg(10)= posInf
      reg(11)= posInf
      reg(12)= posInf

      analyresult=pi**3
      
      do ncall=100, 400, 10
         call vegas(reg,ndim,fcn,0,ncall,itmx,nprn,result,sdr,chi)
         !reg: the integral limit
         !ndim: dimension
         !fcn: grand function
         !0: the initial value
         !ncall: number of throwing
         !itmax: number of interation
         !nprn: output requirement
         !result: the value of integral
         !sdr: standard deviation
         error=(result-analyresult)/analyresult
         write (1,*) ncall, result,error
      end do
      stop
      end program vegasIntegral

      INCLUDE 'vegas.f'

      real*8 function fcn(p)
       implicit none
       real*8 p(6),pi
       pi = 4.d0*datan(1.d0)
       fcn = exp(-p(1)**2-p(2)**2-p(3)**2-p(4)**2-p(5)**2-p(6)**2)
       return
      end

