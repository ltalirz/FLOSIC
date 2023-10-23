C THIS IS A COLLECTION OF USEFUL SUBROUTINES THAT WERE
C DEVELOPED AT THE UNIVERSITY OF WISCONSIN BY CC LIN, 
C RA HEATON, JG HARRISON, K MEDNICK, S. CHUNG AND 
C OTHERS PRIOR TO 1981.
      SUBROUTINE FMTCAL(XX,S0,S2,S4,S6,S8)
C DEVELOPED BY S
      IMPLICIT  REAL*8 (A-H,O-Z)
      LOGICAL ONCE
      DIMENSION F(10,161)
      DIMENSION G(20)
      SAVE
      DATA ONCE/.FALSE./
      DATA ROOTPI        /1.77245385090551603D+00/
      DATA C1,C2,C3,C4,C5/1.33333333333333333D+00,
     &                    0.53333333333333333D+00,
     &                    0.15238095238095238D+00,
     &                    3.28125000000000000D+00,
     &                    0.93750000000000000D+00/
      DATA F(1,1)/0.0D0/
C
C        DIFFERENCE BETWEEN DOUBLE $ SINGLE  ARE
C        REAL*8 CARD  $  FORMAT 55  $  COUPOT WORD SIZE
C
      IF(ONCE)GO TO 907
      ONCE=.TRUE.
      EPS=1.0D-12
      DO 110 I=1,10
  110 F(I,1)=1.0D0/(2*I-1)
      T=0.0D0
      DO 200 L=2,161
      T=0.1D0*(L-1)
      X=2.0D0*T
      DEN=39.0D0
      TERM=1.0D0/DEN
      SUM=TERM
      DO 130 I=2,100
      DEN=DEN+2.0D0
      TERM=TERM*X/DEN
      SUM=SUM+TERM
      Q=TERM/SUM
      IF(Q-EPS) 140,140,130
  130 CONTINUE
  140 EX=EXP(-T)
      G(20)=EX*SUM
C
C        USE DOWNWARD RECURSION
C
      DO 150 I=1,19
      K=21-I
      KK=K-1
  150 G(KK)=(X*G(K)+EX)/(2*KK-1)
C
      DO 160 I=1,10
  160 F(I,L)=G(I)
  200 CONTINUE
  907 CONTINUE
C
      T=XX
      IF (T-1.6D+01) 10,10,40
   10 IT=INT(10*T+1.5D0)
      DT=T-0.1D0*(IT-1)
      DO 30 MP1=1,5
      FMT=F(MP1,IT)
      TERM=1.0D0
      DO 20 I=1,5
      TERM=-TERM*DT/I
   20 FMT=FMT+TERM*F(MP1+I,IT)
      GO TO (21,22,23,24,25),MP1
   21 S0=FMT
      GO TO 26
   22 S2=FMT
      GO TO 26
   23 S4=FMT
      GO TO 26
   24 S6=FMT
      GO TO 26
   25 S8=FMT
   26 CONTINUE
   30 CONTINUE
      GO TO 100
C
C CHANGED SIMPLY TO EXP - POREZAG FEB 95
C
   40 IF(T .GT. 80.0D0) GO TO 1090
      POL1=EXP(-T)
      GO TO 1091
1090  POL1=0.0D0
1091  RECT=1.0D0/T
      RT=SQRT(T)
      RTERM=RECT+2.0D0
      PTERM=ROOTPI/RT
      DO 51 MP1=1,5
      GO TO (50,60,70,80,90),MP1
   90 POL2=((RTERM*RECT+C1)*RECT+C2)*RECT+C3
      FMT=C4*RECT*(PTERM/(T*T*T)-POL1*POL2)
      S8=FMT
      GO TO 51
   80 POL2=(RTERM*RECT+C1)*RECT+C2
      FMT=C5*RECT*(PTERM/(T*T)-POL1*POL2)
      S6=FMT
      GO TO 51
   70 POL2=RTERM*RECT+C1
      FMT=0.375D0*RECT*(PTERM/T-POL1*POL2)
      S4=FMT
      GO TO 51
   60 FMT=0.25D0*RECT*(PTERM-POL1*RTERM)
      S2=FMT
      GO TO 51
   50 FMT=0.5D0*RECT*(ROOTPI*RT-POL1)
      S0=FMT
   51 CONTINUE
  100 RETURN
      END
C
      SUBROUTINE KNMXSF(A1,A2,R,T,WK)
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL DYES
      DIMENSION R(3),T(3),X(3),WK(ND,ND)
      SAVE
      DATA DYES/.FALSE./
      X(1)=T(1)-R(1)
      X(2)=T(2)-R(2)
      X(3)=T(3)-R(3)
      A=A1*A2/(A1+A2)
      AS=A*A
      AC=AS*A
      PI=3.14159265358979324D0
      B=PI/(A1+A2)
      B=SQRT(B)
      B=B**3
      DO 100 I=1,ND
      DO 100 J=1,ND
  100 WK(I,J)=0.0D0
      XS=X(1)**2+X(2)**2+X(3)**2
      IF( XS.LT.1.0D-14 ) GO TO 150
      C=A*XS
      C=EXP(-C)
C THIS IS SS
      WK(1,1)=A*B*C*(3.0D0-2*A*XS)
      DO 105 I=1,3
C THIS IS SPX
      WK(1,I+1)=-AS*B*C*(5.0D0-2*A*XS)*X(I)/A2
C THIS IS PXS
      WK(I+1,1)= AS*B*C*(5.0D0-2*A*XS)*X(I)/A1
      IF(DYES) GO TO 1000
C THIS IS SDXX
      WK(1,I+4)=0.5D0*A*B*C*(3*A2-5*A+2*A*(A-A2)*XS+14*AS*X(I)**2
     &         -4*AC*XS*X(I)**2)/(A2**2)
C THIS IS DXXS
      WK(I+4,1)=0.5D0*A*B*C*(3*A1-5*A+2*A*(A-A1)*XS+14*AS*X(I)**2
     &         -4*AC*XS*X(I)**2)/(A1**2)
 1000 CONTINUE
C THIS IS PXPX
      WK(I+1,I+1)=AS*B*C*(5.0D0-2*A*XS-14*A*X(I)**2+4*AS*XS*X(I)**2) 
     &           /(2*A1*A2)
      IF(DYES) GO TO 1001
C THIS IS DXXPX
      WK(I+4,I+1)=-AS*B*C*(5*A1-21*A+18*AS*X(I)**2+2*A*(3*A-A1)*XS
     &           -4*AC*XS*X(I)**2)*X(I)/(2*A1**2*A2)
C THIS IS PXDXX
      WK(I+1,I+4)= AS*B*C*(5*A2-21*A+18*AS*X(I)**2+2*A*(3*A-A2)*XS
     &           -4*AC*XS*X(I)**2)*X(I)/(2*A1*A2**2)
C THIS IS DXXDXX
      WKW=2*A1*A2-21*AS-(14*A*A1*A2-108*AC)*X(I)**2+6*AC*XS
     &   +(4*AS*A1*A2-24*A**4)*XS*X(I)**2-44*A**4*X(I)**4
     &   +8*A**5*XS*X(I)**4
      WK(I+4,I+4)=-A*B*C*WKW/(4*A1**2*A2**2)
 1001 CONTINUE
  105 CONTINUE
      IF(DYES) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      WK(1,K)=AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A2**2)
C THIS IS DXYS
      WK(K,1)=AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A1**2)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      WK(I+1,J+1)=-AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A1*A2)
      IF(DYES) GO TO 1003
C THIS IS PXDYY
      WK(I+1,J+4)= AS*B*C*X(I)*(5*A2-7*A-2*A*(A2-A)*XS+18*AS*X(J)**2
     &           -4*AC*XS*X(J)**2)/(2*A1*A2**2)
C THIS IS DXXPY
      WK(I+4,J+1)=-AS*B*C*X(J)*(5*A1-7*A-2*A*(A1-A)*XS+18*AS*X(I)**2
     &           -4*AC*XS*X(I)**2)/(2*A1**2*A2)
C THIS IS DXXDYY
      WKW=7*AS-2*A1*A2+2*AS*(7*A2-9*A)*X(I)**2+2*AS*(7*A1-9*A)*X(J)**2
     &   -2*AC*XS+4*AC*(A-A2)*XS*X(I)**2+4*AC*(A-A1)*XS*X(J)**2
     &   +44*A**4*X(I)**2*X(J)**2-8*A**5*XS*X(I)**2*X(J)**2
      WK(I+4,J+4)=A*B*C*WKW/(4*A1**2*A2**2)
 1003 CONTINUE
  115 CONTINUE
      IF(DYES) GO TO 1004
      DO 130 I=1,3
      DO 130 J=1,2
      DO 130 K=2,3
      IF (J.EQ.K) GO TO 130
      IF (I.NE.J.AND.I.NE.K) GO TO 120
      IF (J.EQ.I) GO TO 118
      L=J
      GO TO 119
  118 L=K
  119 CONTINUE
C THIS IS PXDXY
      WK(I+1,J+K+5)=-AC*B*C*(7.0D0-18*A*X(I)**2
     &             -2*A*XS+4*AS*XS*X(I)**2)*X(L)/(2*A1*A2**2)
C THIS IS DXYPX
      WK(J+K+5,I+1)= AC*B*C*(7.0D0-18*A*X(I)**2
     &             -2*A*XS+4*AS*XS*X(I)**2)*X(L)/(2*A1**2*A2)
C THIS IS DXXDXY
      WK(I+4,J+K+5)=AC*B*C*(7*A1-27*A+22*AS*X(I)**2-2*A*(A1-3*A)*XS
     &             -4*AC*XS*X(I)**2)*X(I)*X(L)/(2*A1**2*A2**2)
C THIS IS DXYDXX
      WK(J+K+5,I+4)=AC*B*C*(7*A2-27*A+22*AS*X(I)**2-2*A*(A2-3*A)*XS
     &             -4*AC*XS*X(I)**2)*X(I)*X(L)/(2*A1**2*A2**2)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      WK(I+1,J+K+5)= A**4*B*C*X(I)*X(J)*X(K)*(9.0D0-2*A*XS)/(A1*A2**2) 
C THIS IS DYZPX
      WK(J+K+5,I+1)=-A**4*B*C*X(I)*X(J)*X(K)*(9.0D0-2*A*XS)/(A1**2*A2)
C THIS IS DXXDYZ
      WK(I+4,J+K+5)=AC*B*C*X(J)*X(K)*(7*A1-9*A+22*AS*X(I)**2
     &             +2*A*(A-A1)*XS-4*AC*XS*X(I)**2)/(2*A1**2*A2**2)
C THIS IS DYZDXX
      WK(J+K+5,I+4)=AC*B*C*X(J)*X(K)*(7*A2-9*A+22*AS*X(I)**2
     &             +2*A*(A-A2)*XS-4*AC*XS*X(I)**2)/(2*A1**2*A2**2)
  130 CONTINUE
      DO 146 I=1,2
      DO 146 J=2,3
      IF (I.EQ.J) GO TO 146
      DO 145 K=1,2
      DO 145 L=2,3
      IF (K.EQ.L) GO TO 145
      IF (I.EQ.K.AND.J.EQ.L) GO TO 140
      IF (I.EQ.K) GO TO 131
      IF (I.EQ.L) GO TO 132
      IF (J.EQ.K) GO TO 133
      IF  (J.EQ.L) GO TO 134
  131 M1=I
      M2=J
      M3=L
      GO TO 135
  132 M1=I
      M2=J
      M3=K
      GO TO 135
  133 M1=J
      M2=I
      M3=L
      GO TO 135
  134 M1=J
      M2=I
      M3=K
      GO TO 135
C THIS IS DXYDYZ
  135 CONTINUE
      WK(I+J+5,K+L+5)=-A**4*B*C*(9.0D0-2*A*XS-22*A*X(M1)**2
     &               +4*AS*XS*X(M1)**2)*X(M2)*X(M3)/(2*A1**2*A2**2)
      GO TO 145
C THIS IS DXYDXY
  140 CONTINUE
      WK(I+J+5,K+L+5)= AC*B*C*(7.0D0-2*A*XS+(4*AS*XS-18*A)*(X(I)**2
     &               +X(J)**2)+(44*AS-8*AC*XS)*X(I)**2*X(J)**2)
     &               /(4*A1**2*A2**2)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      GO TO 200
  150 CONTINUE
C THIS SS
      WK(1,1)=3.0D0*A*B
      DO 155 I=1,3
      IF(DYES) GO TO 1005
C THIS IS SDXX
      WK(1,I+4)=0.5D0*A*B*(3.0D0-5*A/A2)/A2
C THIS IS DXXS
      WK(I+4,1)=0.5D0*A*B*(3.0D0-5*A/A1)/A1
 1005 CONTINUE
C THIS IS PXPX
      WK(I+1,I+1)= 5*B*AS/(2*A1*A2)
      IF(DYES) GO TO 1006
C THIS IS DXXDXX
      WK(I+4,I+4)=-B*A*(2.0D0-21*AS/(A1*A2))/(4*A1*A2)
 1006 CONTINUE
  155 CONTINUE
      IF(DYES) GO TO 1007
      DO 165 I=1,3
      DO 165 J=1,3
      IF (I.EQ.J) GO TO 165
C THIS IS DXXDYY
      WK(I+4,J+4)=-B*A*(2.0D0-7*AS/(A1*A2))/(4*A1*A2)
  165 CONTINUE
      DO 175 I=8,10
C THIS IS DXYDXY
  175 WK(I,I)=7*A*AS*B/(4*A1**2*A2**2)
 1007 CONTINUE
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE OVMXSF(A1,A2,R,T,WO)
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL DYES
      DIMENSION WO(ND,ND),R(3),T(3),X(3)
      SAVE
      DATA DYES/.FALSE./
      X(1)=T(1)-R(1)
      X(2)=T(2)-R(2)
      X(3)=T(3)-R(3)
      A=A1*A2/(A1+A2)
      AS=A*A
      PI=3.14159265358979324D0
      B=PI/(A1+A2)
      B=SQRT(B)
      B=B**3
      DO 100 I=1,ND
      DO 100 J=1,ND
  100 WO(I,J)=0.0D0
      XS=X(1)**2+X(2)**2+X(3)**2
      IF( XS.LT.1.0D-14 ) GO TO 150
      C=A*XS
      C=EXP(-C)
C THIS SS
      WO(1,1)=B*C
      DO 105 I=1,3
C THIS IS SPX
      WO(1,I+1)=-A*B*C*X(I)/A2
C THIS IS PXS
      WO(I+1,1)=A*B*C*X(I)/A1
      IF(DYES) GO TO 1000
C THIS IS SDXX
      WO(1,I+4)=(0.5D0*B*C/A2)*((1.0D0-A/A2)+(2*A*A*X(I)*X(I)/A2))
C THIS IS DXXS
      WO(I+4,1)=(0.5D0*B*C/A1)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
 1000 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=A*B*C*(0.5D0-A*X(I)*X(I))/(A1*A2)
      IF(DYES) GO TO 1001
C THIS IS DXXPX
      WO(I+4,I+1)=-A*B*C*X(I)*((1.0D0-3*A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS PXDXX
      WO(I+1,I+4)= A*B*C*X(I)*((1.0D0-3*A/A2)+(2*A*A*X(I)*X(I)/A2))
     &           /(2*A1*A2)
C THIS IS DXXDXX
      WOW=3.0D0+X(I)**2*(2*(A1+A2)-12*A)+4*AS*X(I)**4
      WO(I+4,I+4)=WOW*B*C*AS/(4*A1**2*A2**2)
 1001 CONTINUE
  105 CONTINUE
      IF(DYES) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      WO(1,K)=A*A*B*C*X(I)*X(J)/(A2*A2)
C THIS IS DXYS
      WO(K,1)=A*A*B*C*X(I)*X(J)/(A1*A1)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      WO(I+1,J+1)=-A*A*B*C*X(I)*X(J)/(A1*A2)
      IF(DYES) GO TO 1003
C THIS IS PXDYY
      WO(I+1,J+4)=A*B*C*X(I)*((1.0D0-A/A2)+(2*A*A*X(J)*X(J)/A2))
     &           /(2*A1*A2)
C THIS IS DXXPY
      WO(I+4,J+1)=-A*B*C*X(J)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS DXXDYY
      WOW=1.0D0+2*(A1-A)*X(J)**2+2*(A2-A)*X(I)**2+4*AS*X(I)**2*X(J)**2
      WO(I+4,J+4)=WOW*B*C*AS/(4*A1**2*A2**2)
 1003 CONTINUE
  115 CONTINUE
      IF(DYES) GO TO 1004
      DO 130 I=1,3
      DO 130 J=1,2
      DO 130 K=2,3
      IF (J.EQ.K) GO TO 130
      IF (I.NE.J.AND.I.NE.K) GO TO 120
      IF (J.EQ.I) GO TO 118
      L=J
      GO TO 119
  118 L=K
  119 CONTINUE
C THIS IS PXDXY
      WO(I+1,J+K+5)=-A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A2*A2)
C THIS IS DXXDXY
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPX
      WO(J+K+5,I+1)=A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A1*A2)
C THIS IS DXYDXX
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      WO(I+1,J+K+5)=(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A2*A2)
C THIS IS DXXDYZ
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPZ
      WO(J+K+5,I+1)=-(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A1*A2)
C THIS IS DXYDZZ
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
  130 CONTINUE
      DO 146 I=1,2
      DO 146 J=2,3
      IF (I.EQ.J) GO TO 146
      DO 145 K=1,2
      DO 145 L=2,3
      IF (K.EQ.L) GO TO 145
      IF (I.EQ.K.AND.J.EQ.L) GO TO 140
      IF (I.EQ.K) GO TO 131
      IF (I.EQ.L) GO TO 132
      IF (J.EQ.K) GO TO 133
      IF  (J.EQ.L) GO TO 134
  131 M1=I
      M2=J
      M3=L
      GO TO 135
  132 M1=I
      M2=J
      M3=K
      GO TO 135
  133 M1=J
      M2=I
      M3=L
      GO TO 135
  134 M1=J
      M2=I
      M3=K
      GO TO 135
C THIS IS DXYDYZ
  135 WO(I+J+5,K+L+5)=-(A**3)*B*C*X(M2)*X(M3)
     &               *(1.0D0-2*A*X(M1)*X(M1))/(2*A1*A1*A2*A2)
      GO TO 145
C THIS IS DXYDXY
  140 WO(I+J+5,K+L+5)=A*A*B*C*(1.0D0-2*A*X(I)*X(I)-2*A*X(J)*X(J)
     &               +4*A*A*X(I)*X(I)*X(J)*X(J))/(4*A1*A1*A2*A2)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      GO TO 200
  150 CONTINUE
C THIS SS
      WO(1,1)=B
      DO 155 I=1,3
      IF(DYES) GO TO 1005
C THIS IS SDXX
      WO(1,I+4)=B*(1.0D0-A/A2)/(2*A2)
C THIS IS DXXS
      WO(I+4,1)=B*(1.0D0-A/A1)/(2*A1)
 1005 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=0.5D0*A*B/(A1*A2)
      IF(DYES) GO TO 1006
C THIS IS DXXDXX
      WO(I+4,I+4)=3*A**2*B/(4*A1**2*A2**2)
 1006 CONTINUE
  155 CONTINUE
      IF(DYES) GO TO 1007
      DO 165 I=1,3
      DO 165 J=1,3
      IF (I.EQ.J) GO TO 165
C THIS IS DXXDYY
      WO(I+4,J+4)=A*A*B/(4*A1**2*A2**2)
C THIS IS DXYDXY
      WO(I+J+5,I+J+5)=WO(I+4,J+4)
  165 CONTINUE
 1007 CONTINUE
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE THCNDR(A1,A2,A3,B,C,DR)
      PARAMETER (ND=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL NOD
      DIMENSION F(3),AE(3),BE(3),DR(ND,ND),B(3),C(3)
      SAVE
      DATA PI/3.141592653589793D0/
      DATA NOD/.FALSE./
      DO 5 K=1,ND
      DO 5 J=1,ND
    5 DR(J,K)=0.0D0
      AT=A1+A2+A3
      EX=(A2*B(1)+A3*C(1))/AT
      EY=(A2*B(2)+A3*C(2))/AT
      EZ=(A2*B(3)+A3*C(3))/AT
      E2=EX**2+EY**2+EZ**2
      B2=B(1)**2+B(2)**2+B(3)**2
      C2=C(1)**2+C(2)**2+C(3)**2
      CCC=AT*E2-A2*B2-A3*C2
      COEF=2.0D0*PI*EXP(CCC)/AT
      F(1)=(A2*B(1)-(A1+A2)*C(1))/AT
      F(2)=(A2*B(2)-(A1+A2)*C(2))/AT
      F(3)=(A2*B(3)-(A1+A2)*C(3))/AT
      F2=F(1)**2+F(2)**2+F(3)**2
      AE(1)=EX
      AE(2)=EY
      AE(3)=EZ
      BE(1)=EX-B(1)
      BE(2)=EY-B(2)
      BE(3)=EZ-B(3)
      XX=AT*F2
      CALL FMTCAL(XX,S0,S2,S4,S6,S8)
C THIS SS
      DR(1,1)=COEF*S0
      DO 105 I=1,3
C THIS IS SPX
      DR(1,I+1)=COEF*(BE(I)*S0-F(I)*S2)
C THIS IS PXS
      DR(I+1,1)=COEF*(AE(I)*S0-F(I)*S2)
      IF(NOD) GO TO 1000
C THIS IS SDXX
      DR(1,I+4)=COEF*((BE(I)**2+0.5D0/AT)*S0-(2*BE(I)*F(I)
     &         +0.5D0/AT)*S2+F(I)**2*S4)
C THIS IS DXXS
      DR(I+4,1)=COEF*((AE(I)**2+0.5D0/AT)*S0-(2*AE(I)*F(I)
     &         +0.5D0/AT)*S2+F(I)**2*S4)
 1000 CONTINUE
C THIS IS PXPX
      DR(I+1,I+1)=COEF*((AE(I)*BE(I)+0.5D0/AT)*S0-(AE(I)*F(I)
     &           +BE(I)*F(I)+0.5D0/AT)*S2+F(I)**2*S4)
      IF(NOD) GO TO 1001
C THIS IS DXXPX
      DP=(AE(I)**2*BE(I)+0.5D0*(2*AE(I)+BE(I))/AT)*S0
      DP=DP-(AE(I)**2*F(I)+0.5D0*(3*F(I)+2*AE(I)+BE(I))/AT
     &  +2*AE(I)*BE(I)*F(I))*S2
      DP=DP+(F(I)**2*(2*AE(I)+BE(I))+1.5D0*F(I)/AT)*S4-F(I)**3*S6
      DR(I+4,I+1)=COEF*DP
C THIS IS PXDXX
      PD=(AE(I)*BE(I)**2+0.5D0*(AE(I)+2*BE(I))/AT)*S0
      PD=PD-(BE(I)**2*F(I)+0.5D0*(3*F(I)+2*BE(I)+AE(I))/AT
     &  +2*AE(I)*BE(I)*F(I))*S2
      PD=PD+(F(I)**2*(AE(I)+2*BE(I))+1.5D0*F(I)/AT)*S4-F(I)**3*S6
      DR(I+1,I+4)=COEF*PD
C THIS IS DXXDXX
      DD=(AE(I)**2*BE(I)**2+0.5D0*(AE(I)**2+BE(I)**2+4*AE(I)*BE(I))/AT
     &  +0.75D0/AT**2)*S0
      DD=DD-((2*AE(I)*BE(I)+3.0D0/AT)*(AE(I)+BE(I))*F(I)
     &  +0.5D0*(AE(I)**2+BE(I)**2+4*BE(I)*AE(I))/AT+1.5D0/AT**2)*S2
      DD=DD+((AE(I)**2+BE(I)**2+4*AE(I)*BE(I))*F(I)**2
     &  +3*F(I)*(AE(I)+BE(I)+F(I))/AT+0.75D0/AT**2)*S4
      DD=DD-(2*F(I)**3*(AE(I)+BE(I))+3*F(I)**2/AT)*S6+F(I)**4*S8
      DR(I+4,I+4)=COEF*DD
 1001 CONTINUE
  105 CONTINUE
      IF(NOD) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      DR(1,K)=COEF*(BE(I)*BE(J)*S0-(BE(J)*F(I)+BE(I)*F(J))*S2+F(I)
     &       *F(J)*S4)
C THIS IS DXYS
      DR(K,1)=COEF*(AE(I)*AE(J)*S0-(AE(J)*F(I)+AE(I)*F(J))*S2+F(I)
     &       *F(J)*S4)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      DR(I+1,J+1)=COEF*(AE(I)*BE(J)*S0-(AE(I)*F(J)+BE(J)*F(I))*S2
     &           +F(I)*F(J)*S4)
      IF(NOD) GO TO 1003
C THIS IS PXDYY
      PD=(BE(J)**2+0.5D0/AT)*AE(I)*S0-F(J)**2*F(I)*S6
      PD=PD-(2*BE(J)*F(J)*AE(I)+0.5D0*AE(I)/AT+BE(J)**2*F(I)
     &  +0.5D0*F(I)/AT)*S2
      PD=PD+(F(J)**2*AE(I)+2*BE(J)*F(J)*F(I)+0.5D0*F(I)/AT)*S4
      DR(I+1,J+4)=COEF*PD
C THIS IS DXXPY
      DP=(AE(I)**2+0.5D0/AT)*BE(J)*S0-F(I)**2*F(J)*S6
      DP=DP-(2*AE(I)*F(I)*BE(J)+0.5D0*BE(J)/AT+AE(I)**2*F(J)
     &  +0.5D0*F(J)/AT)*S2
      DP=DP+(F(I)**2*BE(J)+2*AE(I)*F(I)*F(J)+0.5D0*F(J)/AT)*S4
      DR(I+4,J+1)=COEF*DP
C THIS IS DXXDYY
      DD=(AE(I)**2+0.5D0/AT)*(BE(J)**2+0.5D0/AT)*S0
      DD=DD-(2*AE(I)*BE(J)*(BE(J)*F(I)+AE(I)*F(J))+0.5D0*(AE(I)**2
     &  +BE(J)**2)/AT+(AE(I)*F(I)+BE(J)*F(J))/AT+0.5D0/AT**2)*S2
      DD=DD+(F(I)**2*BE(J)**2+F(J)**2*AE(I)**2+4*AE(I)*BE(J)*F(I)
     &  *F(J)+(AE(I)*F(I)+BE(J)*F(J))/AT+0.5D0*(F(I)**2+F(J)**2)/AT
     &  +0.25D0/AT**2)*S4
      DD=DD-(2*F(I)*F(J)*(BE(J)*F(I)+AE(I)*F(J))
     &  +0.5D0*(F(I)**2+F(J)**2)/AT)*S6
      DD=DD+F(I)**2*F(J)**2*S8
      DR(I+4,J+4)=COEF*DD
 1003 CONTINUE
  115 CONTINUE
      IF(NOD) GO TO 1004
      DO 130 I=1,3
      DO 130 J=1,2
      DO 130 K=2,3
      IF (J.EQ.K) GO TO 130
      IF (I.NE.J.AND.I.NE.K) GO TO 120
      IF (J.EQ.I) GO TO 118
      L=J
      GO TO 119
  118 L=K
  119 CONTINUE
C THIS IS PXDXY
      PD=(BE(L)*S0-F(L)*S2)*(AE(I)*BE(I)+0.5D0/AT)
      PD=PD+(F(L)*S4-BE(L)*S2)*(AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)
      PD=PD+(BE(L)*S4-F(L)*S6)*F(I)**2
      DR(I+1,J+K+5)=COEF*PD
C THIS IS DXXDXY
      DD=(BE(L)*S0-F(L)*S2)*(AE(I)**2*BE(I)+0.5D0*(2*AE(I)+BE(I))/AT)
      DD=DD+(F(L)*S4-BE(L)*S2)*(AE(I)**2*F(I)+0.5D0*(3*F(I)
     &  +2*AE(I)+BE(I))/AT+2*AE(I)*BE(I)*F(I))
      DD=DD+(BE(L)*S4-F(L)*S6)*(F(I)**2*(2*AE(I)+BE(I))+1.5D0*F(I)/AT)
      DD=DD+(F(L)*S8-BE(L)*S6)*F(I)**3
      DR(I+4,J+K+5)=COEF*DD
C THIS IS DXYPX
      DP=(AE(L)*S0-F(L)*S2)*(AE(I)*BE(I)+0.5D0/AT)
      DP=DP+(F(L)*S4-AE(L)*S2)*(AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)
      DP=DP+(AE(L)*S4-F(L)*S6)*F(I)**2
      DR(J+K+5,I+1)=COEF*DP
C THIS IS DXYDXX
      DD=(AE(L)*S0-F(L)*S2)*(AE(I)*BE(I)**2+0.5D0*(AE(I)+2*BE(I))/AT)
      DD=DD+(F(L)*S4-AE(L)*S2)*(BE(I)**2*F(I)+0.5D0*(3*F(I)
     &  +2*BE(I)+AE(I))/AT+2*AE(I)*BE(I)*F(I))
      DD=DD+(AE(L)*S4-F(L)*S6)*(F(I)**2*(AE(I)+2*BE(I))+1.5D0*F(I)/AT)
      DD=DD+(F(L)*S8-AE(L)*S6)*F(I)**3
      DR(J+K+5,I+4)=COEF*DD
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      PD=AE(I)*BE(J)*BE(K)*S0-(AE(I)*(BE(J)*F(K)+BE(K)*F(J))
     &  +F(I)*BE(J)*BE(K))*S2
      PD=PD+(AE(I)*F(J)*F(K)+BE(K)*F(I)*F(J)+BE(J)*F(I)*F(K))*S4
     &  -F(I)*F(J)*F(K)*S6
      DR(I+1,J+K+5)=COEF*PD
C THIS IS DXXDYZ
      DD=(BE(K)*S0-F(K)*S2)*((AE(I)**2+0.5D0/AT)*BE(J))
      DD=DD+(F(K)*S4-BE(K)*S2)*(2*AE(I)*F(I)*BE(J)+0.5D0*BE(J)/AT
     &  +AE(I)**2*F(J)+0.5D0*F(J)/AT)
      DD=DD+(BE(K)*S4-F(K)*S6)*(F(I)**2*BE(J)+2*AE(I)*F(I)*F(J)
     &  +0.5D0*F(J)/AT)
      DD=DD+(F(K)*S8-BE(K)*S6)*(F(I)**2*F(J))
      DR(I+4,J+K+5)=COEF*DD
C THIS IS DYZPX
      DP=AE(J)*AE(K)*BE(I)*S0-(F(J)*AE(K)*BE(I)+AE(J)*F(K)*BE(I)
     &  +AE(J)*AE(K)*F(I))*S2
      DP=DP+(F(J)*F(K)*BE(I)+F(J)*AE(K)*F(I)+AE(J)*F(K)*F(I))*S4
      DP=DP-F(J)*F(K)*F(I)*S6
      DR(J+K+5,I+1)=COEF*DP
C THIS IS DYZDXX
      DD=(AE(K)*S0-F(K)*S2)*(BE(I)**2+0.5D0/AT)*AE(J)
      DD=DD+(F(K)*S4-AE(K)*S2)*(2*BE(I)*F(I)*AE(J)+0.5D0*AE(J)/AT
     &  +BE(I)**2*F(J)+0.5D0*F(J)/AT)
      DD=DD+(AE(K)*S4-F(K)*S6)*(F(I)**2*AE(J)+2*BE(I)*F(I)*F(J)
     &  +0.5D0*F(J)/AT)
      DD=DD+(F(K)*S8-AE(K)*S6)*F(I)**2*F(J)
      DR(J+K+5,I+4)=COEF*DD
  130 CONTINUE
      DO 146 I=1,2
      DO 146 J=2,3
      IF (I.EQ.J) GO TO 146
      DO 145 K=1,2
      DO 145 L=2,3
      IF (K.EQ.L) GO TO 145
      IF (I.EQ.K.AND.J.EQ.L) GO TO 140
      IF (I.EQ.K) GO TO 131
      IF (I.EQ.L) GO TO 132
      IF (J.EQ.K) GO TO 133
      IF  (J.EQ.L) GO TO 134
  131 M1=I
      M2=J
      M3=L
      GO TO 135
  132 M1=I
      M2=J
      M3=K
      GO TO 135
  133 M1=J
      M2=I
      M3=L
      GO TO 135
  134 M1=J
      M2=I
      M3=K
      GO TO 135
  135 CONTINUE
C THIS IS DXYDYZ
      DD=(BE(M3)*S0-F(M3)*S2)*(AE(M1)*BE(M1)+0.5D0/AT)*AE(M2)
      DD=DD+(F(M3)*S4-BE(M3)*S2)*((AE(M1)*F(M1)+BE(M1)*F(M1)
     &  +0.5D0/AT)*AE(M2)+(AE(M1)*BE(M1)+0.5D0/AT)*F(M2))
      DD=DD+(BE(M3)*S4-F(M3)*S6)*(F(M1)**2*AE(M2)+(AE(M1)*F(M1)
     &  +BE(M1)*F(M1)+0.5D0/AT)*F(M2))
      DD=DD+(F(M3)*S8-BE(M3)*S6)*F(M2)*F(M1)**2
      DR(I+J+5,K+L+5)=COEF*DD
      GO TO 145
  140 CONTINUE
C THIS IS DXYDXY
      DD=(AE(I)*BE(I)+0.5D0/AT)*(AE(J)*BE(J)+0.5D0/AT)*S0
      DD=DD-((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*(AE(J)*BE(J)+0.5D0/AT)
     &  +(AE(J)*F(J)+BE(J)*F(J)+0.5D0/AT)*(AE(I)*BE(I)+0.5D0/AT))*S2
      DD=DD+((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*(AE(J)*F(J)+BE(J)*F(J)
     &  +0.5D0/AT)+(AE(I)*BE(I)+0.5D0/AT)*F(J)**2+(AE(J)*BE(J)
     &  +0.5D0/AT)*F(I)**2)*S4
      DD=DD-((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*F(J)**2
     &  +(AE(J)*F(J)+BE(J)*F(J)+0.5D0/AT)*F(I)**2)*S6
      DD=DD+F(I)**2*F(J)**2*S8
      DR(I+J+5,K+L+5)=COEF*DD
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      RETURN
      END
C
      SUBROUTINE THCNOV(A1,A2,A3,B,C,OV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL NOD
      DIMENSION AE(3),BE(3),OV(ND,ND),B(3),C(3)
      SAVE
      DATA NOD/.FALSE./
      DATA PI/3.141592653589793D0/
      DO 5 K=1,ND
      DO 5 J=1,ND
    5 OV(J,K)=0.0D0
      AT=A1+A2+A3
      EX=(A2*B(1)+A3*C(1))/AT
      EY=(A2*B(2)+A3*C(2))/AT
      EZ=(A2*B(3)+A3*C(3))/AT
      E2=EX**2+EY**2+EZ**2
      B2=B(1)**2+B(2)**2+B(3)**2
      C2=C(1)**2+C(2)**2+C(3)**2
      CCC=AT*E2-A2*B2-A3*C2
      COEF=SQRT(PI/AT)**3*EXP(CCC)
      AE(1)=EX
      AE(2)=EY
      AE(3)=EZ
      BE(1)=(EX-B(1))
      BE(2)=(EY-B(2))
      BE(3)=(EZ-B(3))
C THIS SS
      OV(1,1)=COEF
      DO 105 I=1,3
C THIS IS SPX
      OV(1,I+1)=COEF*BE(I)
C THIS IS PXS
      OV(I+1,1)=COEF*AE(I)
      IF(NOD) GO TO 1000
C THIS IS SDXX
      OV(1,I+4)=COEF*(BE(I)**2+0.5D0/AT)
C THIS IS DXXS
      OV(I+4,1)=COEF*(AE(I)**2+0.5D0/AT)
 1000 CONTINUE
C THIS IS PXPX
      OV(I+1,I+1)=COEF*(AE(I)*BE(I)+0.5D0/AT)
      IF(NOD) GO TO 1001
C THIS IS DXXPX
      OV(I+4,I+1)=COEF*(BE(I)*AE(I)**2+0.5D0*(2.0D0*AE(I)+BE(I))/AT)
C THIS IS PXDXX
      OV(I+1,I+4)=COEF*(AE(I)*BE(I)**2+0.5D0*(2.0D0*BE(I)+AE(I))/AT)
C THIS IS DXXDXX
      OV(I+4,I+4)=COEF*(AE(I)**2*BE(I)**2+0.5D0*(AE(I)**2+BE(I)**2) 
     &           /AT+2*BE(I)*AE(I)/AT+0.75D0/AT**2)
 1001 CONTINUE
  105 CONTINUE
      IF(NOD) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      OV(1,K)=COEF*BE(I)*BE(J)
C THIS IS DXYS
      OV(K,1)=COEF*AE(I)*AE(J)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      OV(I+1,J+1)=COEF*AE(I)*BE(J)
      IF(NOD) GO TO 1003
C THIS IS PXDYY
      OV(I+1,J+4)=COEF*AE(I)*(BE(J)**2+0.5D0/AT)
C THIS IS DXXPY
      OV(I+4,J+1)=COEF*BE(J)*(AE(I)**2+0.5D0/AT)
C THIS IS DXXDYY
      OV(I+4,J+4)=COEF*(AE(I)**2+0.5D0/AT)*(BE(J)**2+0.5D0/AT)
 1003 CONTINUE
  115 CONTINUE
      IF(NOD) GO TO 1004
      DO 130 I=1,3
      DO 130 J=1,2
      DO 130 K=2,3
      IF (J.EQ.K) GO TO 130
      IF (I.NE.J.AND.I.NE.K) GO TO 120
      IF (J.EQ.I) GO TO 118
      L=J
      GO TO 119
  118 L=K
  119 CONTINUE
C THIS IS PXDXY
      OV(I+1,J+K+5)=COEF*BE(L)*(AE(I)*BE(I)+0.5D0/AT)
C THIS IS DXXDXY
      OV(I+4,J+K+5)=COEF*BE(L)*(AE(I)**2*BE(I)
     &             +0.5D0*(2*AE(I)+BE(I))/AT)
C THIS IS DXYPX
      OV(J+K+5,I+1)=COEF*AE(L)*(AE(I)*BE(I)+0.5D0/AT)
C THIS IS DXYDXX
      OV(J+K+5,I+4)=COEF*AE(L)*(AE(I)*BE(I)**2
     &             +0.5D0*(2*BE(I)+AE(I))/AT)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      OV(I+1,J+K+5)=COEF*BE(J)*BE(K)*AE(I)
C THIS IS DXXDYZ
      OV(I+4,J+K+5)=COEF*BE(J)*BE(K)*(AE(I)**2+0.5D0/AT)
C THIS IS DYZPX
      OV(J+K+5,I+1)=COEF*AE(J)*AE(K)*BE(I)
C THIS IS DYZDXX
      OV(J+K+5,I+4)=COEF*AE(J)*AE(K)*(BE(I)**2+0.5D0/AT)
  130 CONTINUE
      DO 146 I=1,2
      DO 146 J=2,3
      IF (I.EQ.J) GO TO 146
      DO 145 K=1,2
      DO 145 L=2,3
      IF (K.EQ.L) GO TO 145
      IF (I.EQ.K.AND.J.EQ.L) GO TO 140
      IF (I.EQ.K) GO TO 131
      IF (I.EQ.L) GO TO 132
      IF (J.EQ.K) GO TO 133
      IF  (J.EQ.L) GO TO 134
  131 M1=I
      M2=J
      M3=L
      GO TO 135
  132 M1=I
      M2=J
      M3=K
      GO TO 135
  133 M1=J
      M2=I
      M3=L
      GO TO 135
  134 M1=J
      M2=I
      M3=K
      GO TO 135
  135 CONTINUE
C THIS IS DXYDYZ
      OV(I+J+5,K+L+5)=COEF*AE(M2)*BE(M3)*(AE(M1)*BE(M1)+0.5D0/AT)
      GO TO 145
  140 CONTINUE
C THIS IS DXYDXY
      OV(I+J+5,K+L+5)=COEF*(AE(I)*BE(I)+0.5D0/AT)
     &               *(AE(J)*BE(J)+0.5D0/AT)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      RETURN
      END
C
