program main
IMPLICIT NONE
integer,parameter::L=12,N=4096,L0=924,m=50 !# of particles,# of all possible states, # of Sz=0 states
integer::i,j
double precision::g=1.0d0
double precision,dimension(L0)::Hphi,phi,phiprev,phinew 
double precision::phiHphi,phi2,phiprev2
integer,dimension(L0)::relable

!the parameters in the tridiagonal form of H in phi basis
double precision,dimension(0:m)::a,aa
double precision,dimension(0:m)::b,bb
double precision,dimension(m+1,m+1)::z
integer::flag
z=0d0


call basis(L,N,L0,relable)

phiprev=0d0
do i=1,L0
   if (i==1.or.i==2) then
      phi(i)=1d0/dsqrt(2d0)
   else 
      phi(i)=0d0
end if
end do


a=0d0
b=0d0
print *,'g=',g
phiprev2=0d0
do j=0,m
!loop of phi_j to get a_j, b_j until |phi_j|=0
  
 phiHphi=0d0
   phi2=0d0
   call Hamiltonian(relable,L0,L,g,phi,Hphi)
!   print *,'j=',j
!   print *,'after call,phi=',phi,'Hphi=',Hphi
   do i=1,L0!loop inside phi to calculate norm 
      !call Hamiltonian(relable,L0,L,g,phi,Hphi)
      phiHphi=phiHphi+phi(i)*Hphi(i)
!   print *,'i=',i,phi(i)
   phi2=phi2+phi(i)*phi(i)
!   print *,'i=',i,phi(i),Hphi(i)
!   print *,'phiHphi=',phiHphi,'phi2=',phi2

   !IF (phi2==0d0) EXIT 
   end do
!   print *,'j=',j,'phiHphi=',phiHphi,'phi2=',phi2
   a(j)=phiHphi/phi2
   if (j==0) then
      b(j)=0d0!b_0=0
   else
      b(j)=phi2/phiprev2
   end if
!print *,'j=',j,'phi2=',phi2,'phiprev2=',phiprev2
   phiprev2=phi2
   phinew=Hphi-a(j)*phi-b(j)*phiprev
!   print *,'phinew=',phinew
!   print *,'phi=',phi
   phiprev=phi
   phi=phinew
!   print *,'j=',j
!   print *,'phi=',phi
!   print *,'phinew=',phinew
!   j=j+1

if (j.NE.0) then
   aa=a
   bb=dsqrt(b)
   call imtql2(m+1,m+1,aa,bb,z,flag)
print *,'j=',j,'flag=',flag,'E0/JL=',aa(0)/dble(L),'E1/JL=',aa(1)/dble(L),'delta=',(aa(1)-aa(0))/dble(L)
!print *,aa(4)
!print *, 'aa=',aa
end if
end do
!print *,'g=',g,'E0/JL=',aa(0)/dble(L)


!print *,'a=',a
!print *,'b=',b
pause
end program

INCLUDE 'imtql2.f' 
INCLUDE 'pythag.f'
!go through all possible states to find out Sz=0 states, and relable them



subroutine basis(L,N,L0,relable)
integer::N,L,L0
integer,dimension(L0)::relable
integer,dimension(L)::state

integer::i
integer::sc,rsc! state code and relabel state code
integer::Sz !total Sz
External::dec2bin
!integer:: r
rsc=1 
do sc=0,N-1 
   call dec2bin(sc,L,state) 
   Sz=0 
   do i=1,L
      Sz=Sz+state(i)
   end do !loop of calculating Sz
   if (Sz==L/2) then
      relable(rsc)=sc !rsc is the index of relable array;sc is the content of relable array
      rsc=rsc+1
   end if
end do !loop of going thtough all states

do i=1,L0
!   print *,'relabel=', relable(i)
end do
end subroutine basis

subroutine Hamiltonian(relable,L0,L,g,phi,Hphi)
integer::L0,L
double precision::g
double precision,dimension(L0)::phi,Hphi,Temphi
integer,dimension(L0)::relable

integer::i
integer,dimension(L)::state,flip,pl
integer::sc,rsc! state code and relabel state code
integer::fsc,frsc! flip state code and relabel flip state code
double precision::h0
external::dec2bin,bin2dec,search
!pl is the PBC next
do i=1,L
   if (i==L) then
      pl(i)=1
   else 
      pl(i)=i+1
   end if
end do
!print *,'pl=',pl

Hphi=0d0
do rsc=1,L0
   Temphi=0d0
   if(phi(rsc).NE.0) then
     sc=relable(rsc)
!     print *,'rsc=',rsc,'sc=',sc
     call dec2bin(sc,L,state) !decode state
!print *,state     
     h0=0d0
     do i=1,L
        flip=state
 !       print *,'flip=',flip

        if (state(i)==state(pl(i))) then
           h0=h0+1d0
! print *,'i=',i,'h0=',h0
        else
!print *,'i=',i,'h0=',h0

           h0=h0-1d0
           flip(i)=state(pl(i))
           flip(pl(i))=state(i)
!print *,state
!print *,flip          
           call bin2dec(flip,L,fsc)
!print *,fsc
           call search(relable,L0,fsc,frsc)
!print *,frsc
!Hphi(frsc)=Hphi(frsc)+1
           Temphi(frsc)=0.5d0*g
        end if
     end do!loop of different spins in one state
!     print *,'h0=',h0
     Temphi(rsc)=0.25d0*h0
    ! Hphi(rsc)=Hphi(rsc)+h0
!print *,'phi=',phi
!print *,'Temphi=',Temphi
     Hphi=phi(rsc)*Temphi+Hphi
   end if

end do

!print *,phi
!print *,Hphi
end subroutine Hamiltonian


subroutine search(relable,L0,sc,rsc)
IMPLICIT NONE
integer::i,L0,rsc,sc
integer,dimension(L0)::relable
integer::TOP,BOT,MID,FIND

FIND=0
TOP=1
BOT=L0
do  while(TOP<=BOT.AND.FIND==0)
   MID=(TOP+BOT)/2
   If(sc==relable(MID)) then
      FIND=1
      rsc=MID
   else if(sc<relable(MID)) then
      BOT=MID-1
   else
      TOP=MID+1
   end if
end do

if (find==0) then
    print*,'Error!'
end if
end subroutine search


subroutine dec2bin(dec,L,bin) !In:dec,L(the dimension of array bin); Out:bin
IMPLICIT NONE
integer:: i,l,dec,decTemp
integer,dimension (L):: bin !binary presentation of the basis state
decTemp=dec
do i=L,1,-1
   bin(i)=mod(decTemp,2)
   decTemp=decTemp/2
!   print *,'bin=',bin(i)
end do 
end subroutine dec2bin


subroutine bin2dec(bin,L,dec) !In:bin,L(the dimension of array bin); Out:dec
IMPLICIT NONE
integer,dimension(L):: bin
integer::i,L,dec
integer::decTemp
decTemp=0
do i=L,1,-1
   decTemp=decTemp+bin(i)*(2**(L-i))
end do
   dec=decTemp
end subroutine bin2dec

