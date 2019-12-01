program eigenvalue
IMPLICIT NONE
double precision:: a,b,c,d,root
double precision, EXTERNAL:: pade,pA,pB
integer::cond,kl,kr
write (*,*) 'Input the range for root finding A & B:'
read (*,*) A,B
!EXTERNAL BISEC,f 
call BISEC(pade,A,B,1d-7,c,d,cond,kl,kr)
print *, cond
if (cond==1) then
   print *, 'Wrong range!'
else if (cond==2) then
   root=c
   print *, 'root=', root
else if (cond==3) then
   root=(a+b)/2.0
   print *, 'root=',root
else if (cond==4) then
   print *, 'Not continued!'
end if
end program

    SUBROUTINE BISEC(pade,A,B,Delta,C,D,Cond,KL,KR)
    double precision, PARAMETER:: Big=1d30
    INTEGER Cond,K,KL,KR,Max
    double precision A,B,C,D,Delta,YA,YB,YC
    double precision, EXTERNAL:: pade,pa,pb
    K=0
    KL=0
    KR=0
    YA=pade(pa,pb,A)
    YB=pade(pa,pb,B)
    D=B-A
    Cond=0
    Max=1+INT((dLOG(D)-dLOG(Delta))/dLOG(2d0))
    IF (YA*YB.GE.0) THEN
        Cond=1
        RETURN
    ENDIF
    DO WHILE ((Cond.EQ.0).AND.(K.LT.Max))
        C=(A+B)/2d0
        YC=pade(pa,pb,C)
            IF (YC.EQ.0d0) THEN
                A=C
                B=C
                Cond=2
            ELSEIF ((YB*YC).GT.0d0) THEN
                B=C
                YB=YC
                KR=KR+1
            ELSE
                A=C
                YA=YC
                KL=KL+1
            ENDIF
        K=K+1
        WRITE(1,*) k,a,c,b
    ENDDO
    D=B-A
    IF (D.LT.Delta) THEN
        IF (Cond.NE.2) Cond=3
        IF (dABS(YA-YB).GT.Big) Cond=4
    ENDIF
!    PAUSE
    RETURN
   END SUBROUTINE

