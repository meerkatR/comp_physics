
L0=fact(L)/(fact(L/2)*fact(L/2))
Nstate=2**L


ALLOCATE(phi(L0),phiOld(L0),phiOlder(L0))



IF(ALLOCATED(phi)) DEALLOCATE(phi)
IF(ALLOCATED(phiOld)) DEALLOCATE(phiOld)
IF(ALLOCATED(phiOlder)) DEALLOCATE(phiOlder)

****************************
phi=(/0.1d0,0.1d0,0d0,0d0,0.7d0,0.3d0/)

call basis(L,N,L0,relable)
phiprev=0d0

call Hamiltonian(relable,L0,L,g,phi,Hphi)
print *,'main program',Hphi


*********************************************
phiprev=0d0
phi=
do j=0,m!loop of phi_j to get a_j, b_j
   phiHphi=0d0
   phi2=0d0
   do i=1,L0!loop inside phi to calculate norm 
      call Hamiltonian(relable,L0,L,g,phi,Hphi)
      phiHphi=phiHphi+phi(i)*Hphi(i)
      phi2=phi2+phi(i)*phi(i)
   end do
   a(j)=phiHphi/phi2
   if (j=0) then
      b(j)=0d0!b_0=0
   else
      b(j)=phi2/phiprev2
   end if
   phiprev2=phi2
   phinew=Hphi-a(j)*phi-b(j)**2*phiprev
   phiprev=phi
   phi=phinew
end do!loop of phi_j to get a_j, b_j

