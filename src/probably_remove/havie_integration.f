C     ALGORITHM 400 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 13, NO. 10,
C     P. 622.
      DOUBLE PRECISION FUNCTION HRVINT(F, A, B, MAX, ACC, FAC, MFIN)
C HAVIE INTEGRATION WITH AN EXPANDED RUTISHAUSER-
C TYPE SUMMATION PROCEDURE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(17), U(17), TPREV(17), UPREV(17)
      EXTERNAL F
C TEST FOR MAX GREATER THAN 16
         MUX=MAX
         IF(MAX-16) 10, 10, 5
    5    MUX=16
C INITIALIZATION
   10    ENPT=0.5D0*(F(A)+F(B))
         SUMT=0.0D0
         MFIN=1
         N=1
         H=B-A
         SH=H
C BEGIN REPETITIVE LOOP FROM ORDER 1 TO ORDER MAX
   15    T(1)=H*(ENPT+SUMT)
         SUM=0.D0
         NN=N+N
         EN=NN
         EM=SH/EN
C BEGIN RUTISHAUSER EVALUATION OF RECTANGULAR SUMS
C INITIALIZATION
         IF(NN-16) 20, 20, 25
   20    NZ=NN
         GO TO 30
   25    NZ=16
         IF(NN-256) 30, 30, 35
   30    NA=NN
         GO TO 40
   35    NA=256
         IF(NN-4096) 40, 40, 45
   40    NB=NN
         GO TO 50
   45    NB=4096
C DEVELOPMENT OF RECTANGULAR SUMS
   50    DO 70 KC=1,NN,4096
             SUMB=0.0D0
             KK=KC+NB-1
             DO 65 KB=KC,KK,256
                SUMA=0.0D0
                KKK=KB+NA-1
                DO 60 KA=KB,KKK,16
                   SUMZ=0.0D0
                   KFR=KA+NZ-1
                   DO 55 KZ=KA,KFR,2
                       ZKZ=KZ
   55              SUMZ=SUMZ+F(A+ZKZ*EM)
   60           SUMA=SUMZ+SUMA
   65        SUMB=SUMA+SUMB
   70    SUM=SUMB+SUM
C END OF RUTISHAUSER PROCEDURE
         U(1)=H*SUM
         K=1
C BEGIN EXTRAPOLATION LOOP
   75    FAC=DABS(T(K)-U(K))
         IF(T(K)) 80, 85, 80
C TEST FOR RELATIVE ACCURACY
   80    CONTINUE
         IF(FAC-DABS(ACC*T(K))) 90, 90, 100
C TEST FOR ABSOLUTE ACCURACY WHEN T(K)=0
   85    continue
         IF(FAC-DABS(ACC)) 95, 95, 100
   90    continue
         FAC=FAC/DABS(T(K))
C INTEGRAL EVALUATION BEFORE EXIT
   95    continue
         HRVINT=0.5D0*(T(K)+U(K))
         RETURN
  100    continue
         IF(K-MFIN) 105, 115, 115
  105    AK=K+K
         D=2.D0**AK
         DMA=D-1.D0
C BEGIN EXTRAPOLATION
         T(K+1)=(D*T(K)-TPREV(K))/DMA
         TPREV(K)=T(K)
         U(K+1)=(D*U(K)-UPREV(K))/DMA
         UPREV(K)=U(K)
C END EXTRAPOLATION
         K=K+1
         IF(K-MUX) 75, 110, 110
C END EXTRAPOLATION LOOP
  110    FAC=ABS(T(K)-U(K))
         IF(T(K)) 90, 95, 90
C ORDER IS INCREASED BY ONE
  115    H=0.5D0*H
         SUMT=SUMT+SUM
         TPREV(K)=T(K)
         UPREV(K)=U(K)
         MFIN=MFIN+1
         N=NN
         GO TO 15
C RETURN FOR NEXT ORDER EXTRAPOLATION
         END
