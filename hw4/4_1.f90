program quarterPi
integer::N
double precision:: A,A_1,A_2,B,B_1,B_2,pi,Pade,error

A_1=0d0
A_2=1d0
B_1=1d0
B_2=0d0

pi = 4.d0*datan(1.d0)

do N=1, 30

   if (N==1) then 
      A=A_1+A_2
      B=B_1+B_2
   else
      A=(2*N-1)*A_1+(N-1)**2*A_2
      B=(2*N-1)*B_1+(N-1)**2*B_2
   end if

   A_2=A_1
   A_1=A
   B_2=B_1
   B_1=B
   Pade=A/B
   Error=Pade-pi/4.d0
   write(1,*) N,Pade*4.d0, pi, Error  
end do

end program quarterPi 
