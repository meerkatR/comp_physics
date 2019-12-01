	program main 
        call BISEC(F,-1d0,1d0,1d-7,c,d,cond,kl,kr)
	end program



        REAL FUNCTION F(X)
        F = SIN(X)/COS(X)
    RETURN
    END

    SUBROUTINE BISEC(F,A,B,Delta,C,D,Cond,KL,KR)
    PARAMETER(Big=1E5)
    INTEGER Cond,K,KL,KR,Max
    REAL A,B,C,D,Delta,YA,YB,YC
    EXTERNAL F
    K=0
    KL=0
    KR=0
    YA=F(A)
    YB=F(B)
    D=B-A
    Cond=0
    Max=1+INT((ALOG(D)-ALOG(Delta))/ALOG(2.0))
    IF (YA*YB.GE.0) THEN
        Cond=1
        RETURN
    ENDIF
    WHILE (Cond.EQ.0).AND.(K.LT.Max)
        C=(A+B)/2.0
        YC=F(C)
            IF (YC.EQ.0) THEN
                A=C
                B=C
                Cond=2
            ELSEIF ((YB*YC).GT.0) THEN
                B=C
                YB=YC
                KR=KR+1
            ELSE
                A=C
                YA=YC
                KL=KL+1
            ENDIF
        K=K+1
        WRITE(9,1000) K,A,C,B
    REPEAT
    D=B-A
    IF (D.LT.Delta) THEN
        IF (Cond.NE.2) Cond=3
        IF ((ABS(YA).GT.Big).AND.(ABS(YB).GT.Big)) Cond=4
    ENDIF
    PAUSE
    RETURN
