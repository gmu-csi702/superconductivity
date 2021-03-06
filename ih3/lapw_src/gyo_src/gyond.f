      PROGRAM GYORF
C  For SH3 only      
      IMPLICIT REAL*8(A-H,O-Z)              
C     PROGRAM EPINT  (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE7)
C RELATIVISTIC
C        THIS PROGRAM CALCULATES THE ELECTRON-PHONON INTERACTION
C     MATRIX ELEMENT
C
C        INPUTS-
C        EF=ENERGY AT THE FERMI LEVEL
C        NEF=TOTAL DENSITY OF STATES
C       RS=RADIS IN ATOMIC UNITS NEEDED TO CALCULATE THE PHASE SHIFTS
C        RW=RADIUS OF WIGNER-SEITZ SPHERE FOR  N1,THE FREE SCATTERER
C     DENSITY OF STATES
C        N(L)=DECOMPOSED DENSITY OF STATES FOR DIFFERENT L VALUES
C        A= THE LOGARITHMIC DERIVATIVES
C     DERIV. OF LOG.DER.  D(8)
C
      REAL*8 N,N1,NEF,MUST
      COMMON /DAP/ZNORM(4)
      COMMON DELTA,    DDELTA,EF,N1,RW,RS,A,D
      COMMON R,POT,ANTLH,SHIFT,RSSQ,FRPRSQ,P1,
     1ICHG,KMAX,NUMBER,KSTART,
     2IDENO,LPLACE,JDERPR,LMAX2,NATOM,NMAX
      DIMENSION DELTA(4),Y(4),J(4),DDELTA(4),N1(4),N(4),A(1,4,4),D(4,4)
      DIMENSION R(800,4), POT(800,4), ICHG(10,4), KMAX(4), ANTLH(4), SHI
     1FT(4), RS(4), NUMBER(4), KSTART(10,4), JPTLPR(4), RSSQ(4), FRPRSQ(
     24), IDENO(4), JDERPR(4), LPLACE(4), P1(10,2,2)
C     DIMENSION TITLE (10)
      CHARACTER*40 TITLE
      DIMENSION RDEN(4),spdeta(3)
      DIMENSION  EL(3),ALFASQ(3),AM(3)
      DIMENSION omega(3),OMSQA(3)
C
      OPEN (8,FILE='gyo1.out',BLANK = 'ZERO')
      OPEN (9,FILE='dos.in',BLANK = 'ZERO')
C       FILE 9 IS THE DENSITY OF STATES DATA
      OPEN (7,FILE='pot.in',BLANK = 'ZERO')
C       FILE 7 IS THE POTENTIAL DATA
      OPEN(20,FILE='gyo2.out',BLANK='ZERO')
C
      LLATOM=0
      PI=3.141593

C     READ (9,2000)  XCC
 2000 FORMAT (F10.5)
C     WRITE(8,2550) XCC
 2550 FORMAT (1X,4HXCC=,F10.5)
C     READ (9,11)  IATOMS      
  11  FORMAT (I5)
      READ (9,100)TITLE,IATOMS
      WRITE(8,100)TITLE,IATOMS
 100  FORMAT (A40,I5)
C     IATOMS=2
      DO 1000  II=1,IATOMS
C     READ (9,1)  RW
  1   FORMAT  (F10.6)
      READ (9,4) EF,NEF,N(1),N(2),N(3),N(4)
  4   FORMAT (6F10.5)
 5    FORMAT(1H1,50X,'         INPUT'//50X,'EF =',F10.7//50X,'NEF =',F10
     1.4//50X,'R =',F10.7//50X,'RS =',F10.7//33X,'DECOMPOSED DENITY OF S
     2TATES FOR L VALUES 0,1,2,3'//36X,4(2X,F10.4)//48X,'LOGRITHMIC DERI
     3VATIVES'//36X,4(2X,F10.7)//50X,'THEIR DERIVATIVES'//36X,4(2X,F10.5
     4))
C  DOS PER SPIN
      NEF=NEF/2.
      N(1)=N(1)/2.
      N(2)=N(2)/2.
      N(3)=N(3)/2.
      N(4)=N(4)/2.
      XC=1.0
      CALL RADIAL(XC)
      RW=RS(II)
      WRITE(8,98) XC
  98  FORMAT     (10X,6HXC   =,F10.4)
c     READ (9,20) AM(II),theta,MUST,LLATOM
      READ (9,20) AM(II),omega(II),MUST,omegalog
C     REMEMBER THIS PROGRAM IS FOR SH3
  20  FORMAT (4F10.4)
      IF(LLATOM.EQ.0) LLATOM=1
      ALATOM= FLOAT(LLATOM)
      WRITE(8,50) AM(II)
  50  FORMAT (1H ,10X,6HAMASS=,F10.4)
       WRITE(8,60) OMEGA(II)
 60   FORMAT (1H ,10X,7HOMEGA=,F20.4)
C     WRITE(8,70) OMEGA2
  70  FORMAT (1H ,10X,7HOMEGA2=,F10.4)
C     DOS=2.0*NEF/13.6
C     MUST=0.26*DOS/(1.0+DOS)
C     WRITE(8,80) MUST
  80  FORMAT (1H ,10X,7HMUSTAR=,F10.4)
      CALL ONTO
      NA=NATOM
      write (51,*) (A(1,I,NA),I=1,4)
      CALL DSFE
      write (52,*) (A(1,I,NA),I=1,4)
      SEFPI= SQRT(EF)/PI
      NA=NATOM
      DO 490  L=1,4
  490 N1(L)= FLOAT(2*L-1)*SEFPI*RSSQ(NATOM )*(ZNORM(L)**2)*ABS(D(L,NA))
      WRITE(8,5) EF,NEF,RW,RS(NA),(N(I),I=1,4),(A(1,I,NA),I=1,4),(D(I,NA
     1),I=1,4)
      WRITE(8,400)
  400 FORMAT(1X,"ARE THE VALUES OF TH8 FREE SCATTERERS AND THE
     1RATIOS")
      DO 410  K=1,4
      KI=K-1
      RDEN(K)=N(K)/N1(K)
  410 WRITE(8,420)KI,N1(K),RDEN(K)
  420 FORMAT(20X,I10,2F15.6)
      SUM=0.
      DO 2 L=1,3
      term1=(SIN(DELTA(L+1)-DELTA(L)))**2
      term2=(N(L)*N(L+1))/(N1(L)*N1(L+1))
      term3=2.*L*term1*term2*EF/((PI**2)*NEF)
      term3=term3*48.59
      spdeta(l)=term3
      write(25,501) l-1,term1,term2,term3,rw,ef,nef

 501  format(1x,"l,SIN,DOS-PRODUCT,ETA,rs,ef,nef",I5,6F10.5)
c 501  format(1x,I5,6F10.5)
      X=(2*L*((SIN((DELTA(L+1))-(DELTA(L))))**2)*(N(L))*(N(L+1)))/
     1((N1(L))*(N1(L+1)))
      SUM=SUM+X
      EPIME=EF*SUM/(((3.14159)**2)*(NEF**2))
       EPI=NEF*EPIME
C        OUTPUTS-
      WRITE(8,401)
  401 FORMAT(1X,12HPHASE SHIFTS)
      WRITE(8,10) DELTA
      WRITE(8,402)
  402 FORMAT(1X,27HDERIVATIVES OF PHASE SHIFTS)
      WRITE(8,10)DDELTA
 10   FORMAT (4F12.6)
C        EPIME=ELECTRON-PHONON INTERACTION MATRIX ELEMENT
C        EPI=INTERACTION ELEMENT TIMES THE TOTAL DENSITY OF STATES
C
      WRITE(8,3) EPIME,EPI
 3    FORMAT(1H1,20X,'ELECTRON-PHONON INTERACTION MATRIX ELEMENT =',E17.
     110///12X,'PRODUCT OF INTERACTION ELEMENT AND TOTAL DENSITY OF STAT
     2ES =',E17.10)
      IF (II.EQ.1)  ETA=EPI*48.59
      IF (II.EQ.2)  ETA=EPI*48.59
c     IF (II.EQ.3)  ETA=EPI*48.59
c     IF (II.EQ.4)  ETA=EPI*48.59
      WRITE (8,81) ETA
 81   FORMAT (1H ,10X,4HETA=,F10.4)
      CL=2.7353285E+07
      OMSQA(II)=OMEGA(II)**2
c     OMSQA=(theta)**2/2.0
        EL(II)=EPI*CL/(AM(II)*OMSQA(II))
c     IF (II.EQ.2)  ELAMDA=EPI*CL/(AM*OMSQA)
      WRITE(8,30) EL(II)
  30  FORMAT (1H ,10X,6HLAMDA=,F8.4)
 2    CONTINUE
      WRITE(20,82) TITLE,ETA
      WRITE(21,83) TITLE,EF,NEF,ETA
 82   FORMAT(A40,F10.5)
 83   FORMAT(A40,3F10.5)
C     OMG1(II)=OMEGA1
      CV=5.6466434E+05
      OMEGA(II)=OMEGA(II)/CV
         ALFASQ(II)=ETA/(2.0*AM(II)*OMEGA(II))
C     WRITE(8,300)ALFAS
 300  FORMAT (1H ,10X,6HALFAS=,F10.5)
 1000 CONTINUE
c     IF (IATOMS.EQ.1)  ELAMD=EL(1)
C     IF (IATOMS.EQ.1)  ALFS=ALFASQ(1)
C     IF (IATOMS.EQ.1)  ALLEN=ALFASQ(1)*OMG1(1)
c     IF (IATOMS.EQ.2)  ELAMD=EL(1)+EL(2)
c   For BCC
       ELAMD=EL(1)+3.0*EL(2)
C     IF (IATOMS.EQ.2)  ALFS=ALFASQ(1)+ALFASQ(2)
C     IF (IATOMS.EQ.2)  ALLEN=ALFASQ(1)*OMG1(1)+ALFASQ(2)*OMG1(2)
C     OMEGA3=2.0*ALFS/ELAMD
C     OMEKAP=(2.0*ALLEN/ELAMD)**0.5
C     OMELOG=2.0*OMEGA3-OMEKAP
C  Special to SH3 under pressure from Li et al       
c     OMELOG=1000.0
      WRITE(8,72) EL(1),EL(2), ELAMD
      WRITE(8,73) ALFASQ(1),ALFASQ(2)
  72  FORMAT (1H ,10X,8HLAMDTOT=,3F10.4)
  73  FORMAT (1H ,10X,7HALFASQ=,2F10.4)
C     WRITE(8,71) OMEGA3
  71  FORMAT (1H ,10X,7HOMEGA3=,F11.4) 
C     WRITE(8,99) OMEKAP
 99   FORMAT (1H ,10X,8HOMEKAP =,F10.4)
C     WRITE(8,999) OMELOG
 999  FORMAT (1H ,10X,8HOMELOG =,F10.4)
      RATIO=-1.04*(1.+ELAMD )/(ELAMD -MUST*(1.+0.62*ELAMD ))
      TC=OMEGALOG*EXP(RATIO)/1.20
C     TC=OMEGA3*EXP(RATIO)/1.20
c     OMEGA=theta/SQRT(2.0)
c     TC=theta*EXP(RATIO)/1.45
c     TC=THETA*EXP(RATIO)/1.20
      TCBCS=1.13*THETA*EXP(-1./elamd)
      WRITE (8,40) TC
c     WRITE (8,42) theta
      WRITE (8,43) must
c     WRITE (8,41) TCBCS
  40  FORMAT (1H ,10X,3HTC=,F8.4)
  42  FORMAT (1H ,10X,6HDebye=,F8.4)
  43  FORMAT (1H ,10X,6HMusta=,F8.4)
  41  FORMAT (1H ,10X,6HTCBCS=,F8.4)
      CLOSE(8)
      CLOSE(9)
      CLOSE(7)
      CLOSE(20)
      STOP
      END
      SUBROUTINE DSFE
      IMPLICIT REAL*8(A-H,O-Z)              
C     This program is  changed to incorporate the 10 doublings
      REAL*8 N,J,K
      COMMON/DAP/Z
      COMMON DELTA,    DDELTA,EF,N ,RW,RS,AA,DD
      COMMON R,POT,ANTLH,SHIFT,RSSQ,FRPRSQ,P1,
     1ICHG,KMAX,NUMBER,KSTART,
     2IDENO,LPLACE,JDERPR,LMAX2,NATOM,NMAX
      DIMENSION R(800,4), POT(800,4), ICHG(10,4), KMAX(4), ANTLH(4), SHI
     1FT(4), RS(4), NUMBER(4), KSTART(10,4), JPTLPR(4), RSSQ(4), FRPRSQ(
     24), IDENO(4), JDERPR(4), LPLACE(4), P1(10,2,2)
      DIMENSION DELTA(4),Y(4),J(4),DDELTA(4),N(4),     DJ(4),DY(4),DR(4)
     1,D(4),Z(4)
      DIMENSION AA(1,4,4),DD(4,4)
      K=SQRT(EF)
      X=K*RW
      C=COS(X)
      S=SIN(X)
      F=-C
C                BESSEL FUNCTIONS
      J(1)=S/X
      J(2)=(S/X**2)-(C/X)
      J(3)=(((3/X**3)-(1./X))*S)-((3./X**2)*C)
      J(4)=(((15./X**4)-(6./X**2))*S)-(((15./X**3)-(1./X))*C)
C DERIVATIVES OF BESSEL FUNCTIONS
      DJ(1)=-S/X**2+C/X
      DJ(2)=(-2.*S)/X**3+C/X**2+C/X**2+S/X
      DJ(3)=(-9.*S)/X**4+(3.*C)/X**3+S/X**2-C/X+(6.*C)/X**3+(3.*S)/X**2
      DJ(4)=-(60.*S)/X**5+(15.*C)/X**4+(12.*S)/X**3-(6.*C)/X**2+(45.*C)/
     .X**4+(15.*S)/X**3-(1./X**2)*C-(1./X)*S
C               NEUMANN FUNCTIONS
      Y(1)=F/X
      Y(2)=(F/X**2)-(S/X)
      Y(3)=(((3/X**3)-(1./X))*F)-((3./X**2)*S)
      Y(4)=(((15./X**4)-(6./X**2))*F)-(((15./X**3)-(1./X))*S)
C DERIVATIVE OF NEUMANN FUNCTIONS
      DY(1)=C/X**2+S/X
      DY(2)=(2.*C)/X**3+S/X**2+S/X**2-C/X
      DY(3)=(9.*C)/X**4+(3.*S)/X**3-(1./X**2)*C-(1./X)*S+(6.*S)/X**3-(3.
     .*C)/X**2
      DY(4)=(60.*C)/X**5+(15.*S)/X**4-(12.*C)/X**3-(6.*S)/X**2+(45.*S)/X
     .**4-(15.*C)/X**3-(1./X**2)*S+(1./X)*C
      DO 27 I=1,4
   27 D(I)=DELTA(I)
      DO 23 L=1,4
      Z(L)=J(L)*COS(D(L))-Y(L)*SIN(D(L))
   23 DR(L)=DJ(L)*COS(D(L))-DY(L)*SIN(D(L))
      PI=3.141593
      DO 3 L=1,4
      A=((2.*(L-1)+1)/PI)
      B=DDELTA(L)
      C=((K/2.)*(RW**3)*(DR(L)**2))
      E=(RW**2/2.)*Z(L)*DR(L)
      F=RW/(2.*K)
      G=(EF*(RW**2))-((L-1)*L)
    3 N(L)=A*(B+C+E+F*G*(Z(L)**2))
      RETURN
      END
      SUBROUTINE INTG (XC)
      IMPLICIT REAL*8(A-H,O-Z)   
      REAL*8 N1           
      COMMON DELTA,    DDELTA,EF,N1,RW,RS ,RSQLND,DLD
      COMMON R,POT,ANTLH,SHIFT,RSSQ,FRPRSQ,P1,
     1ICHG,KMAX,NUMBER,KSTART,
     2IDENO,LPLACE,JDERPR,LMAX2,NATOM,NMAX
      COMMON /FANG/ FREL(800,10),GREL(800,10)
      DIMENSION R(800,4), POT(800,4), ICHG(10,4), KMAX(4), ANTLH(4), SHI
     1FT(4), RS(4), NUMBER(4), KSTART(10,4), JPTLPR(4), RSSQ(4), FRPRSQ(
     24), IDENO(4), JDERPR(4), LPLACE(4), P1(10,2,2)
      DIMENSION         ULUL(1,4,4),PNLSQ(800,4,4),OVER(4,4),DLD(4,4)
      DIMENSION RSQLND(1,4,4)
      DIMENSION DELTA(4),DDELTA(4),N1(4)
      DIMENSION      RLG(20),DX1(800)
      DIMENSION A(7), B(6), RR(800), P(800)
      C=274.07204
C     C=1.E+7
      NN=NATOM
C   IF LDR.NE.1   DOES  VCA  FOR  H
C     READ (9,155) LDR,IDENO(NN),XC
  155 FORMAT (I4,A4,F6.3)
C     WRITE (8,166)  XC
 166  FORMAT (1X, 3HXC=,F6.3)
   20 MAX=KMAX(NN)
      LP=LPLACE(NN)
      M5=LP-3
 15   H=R(LP,NN)-R(LP-1,NN)
      PINTRP=(RS(NN)-R(LP,NN))/H
      HSQ=PINTRP*PINTRP
      PM=PINTRP*(HSQ-1.0)*(HSQ-4.0)*(PINTRP-3.0)
      B(1)=-PM/(120.0*(PINTRP+2.0))
      B(2)=PM/(24.0*(PINTRP+1.0))
      B(3)=-PM/(12.0*PINTRP)
      B(4)=PM/(12.0*(PINTRP-1.0))
      B(5)=-PM/(24.0*(PINTRP-2.0))
      B(6)=PM/(120.0*(PINTRP-3.0))
      RP2= RS(NN)*RS(NN)
      ADX= R(1,NN)
      HOC= POT(1,NN)*ADX/C
      WRITE(60,*) POT(1,NN),ADX,C,HOC
      DO 501  I=1,MAX
      DX1(I)= ADX
      IF(I.GE.32) DX1(I)= 2.*ADX
      IF(I.GE.64) DX1(I)= 4.*ADX
      IF(I.GE.96) DX1(I)= 8.*ADX
      IF(I.GE.128) DX1(I)= 16.*ADX
      IF(I.GE.160) DX1(I)= 32.*ADX
C     Change for 10 doublings 
      IF(I.GE.192) DX1(I)= 64.*ADX
      IF(I.GE.224) DX1(I)=128.*ADX
      IF(I.GE.256) DX1(I)=256.*ADX
      IF(I.GE.288) DX1(I)=512.*ADX 
  501 CONTINUE
      write(61,*) HOC
      CALL RELOGR(EF,MAX,LMAX2 ,B,ADX,DX1,POT(1,NN),RP2,HOC,M5,RLG,1)
      DO 105 L1=1,4
      ML= L1 -1
      GL= 0.0
      GL1= 0.0
      DO 75  LR=1,6
      KR= LR + M5
      GL= GL + GREL(KR,L1)  *B(LR)/R(KR,NN)
   75 CONTINUE
C     GL= GL*GL
   95 ULUL(1,L1,NN)= 1.0
      DO 100 K=1,MAX
  100 PNLSQ(K,L1,NN)= (FREL(K,L1)/GL)**2 + (GREL(K,L1)/GL)**2
      GO TO 160
  150 CONTINUE
 160  CONTINUE
  105 CONTINUE
      DO 125  LL= 1,4
  125 RSQLND(1,LL,NN)= RLG(LL)
  130 CONTINUE
      JDERPR(NN)=1
C     PRINT OUT LOGARITHMIC DERIVATIVES
      WRITE (8,165) NN,IDENO(NN)
      WRITE(8,191)(RSQLND(1,LL,NN),LL=1,4)
  191 FORMAT(4E24.8)
      N=NATOM
      DO 80 L=1,4
      OVER(L,N)=SIMCHG(PNLSQ(1,L,N),R(1,N),N,LPLACE(1),RS(1))
C     DLD IS THE DER. OF LOGDER WITH RESPECT TO ENERGY TIMES RS
      DLD(L,N)=OVER(L,N)/ULUL(1,L,N)
      DLD(L,N)=-DLD(L,N)
      WRITE(8,161)   DLD(L,N)
   80 CONTINUE
C
  161 FORMAT ( E24.8)
  165 FORMAT (1H1,29X,62H FOLLOWING ARE VALUES OF SPHERE RADIUS SQUARED
     1TIMES THE VALUE//27X,67H OF THE LOGARITHMIC DERIVATIVE AT THE SPHE
     2RE RADIUS FOR ATOM NUMBER,I4//49X,24H IDENTIFICATION NUMBER  ,A6//
     3/42X,34HCOLUMNS ARE HEADED BY THE ENERGIES/////)
  170 FORMAT (///5H     ,6E18.8)
  175 FORMAT (1H )
  180 FORMAT (5H     ,6E18.8)
  110 JDERPR(NN)=1
  115 CONTINUE
      RETURN
      END
      SUBROUTINE ONTO
      IMPLICIT REAL*8(A-H,O-Z)              
      REAL*8 K,N1
      COMMON Q,    W,EF,N1,RW,RS,A,DD
      COMMON R,POT,ANTLH,SHIFT,RSSQ,FRPRSQ,P1,
     1ICHG,KMAX,NUMBER,KSTART,
     2IDENO,LPLACE,JDERPR,LMAX2,NATOM,NMAX
      DIMENSION R(800,4), POT(800,4), ICHG(10,4), KMAX(4), ANTLH(4), SHI
     1FT(4), RS(4), NUMBER(4), KSTART(10,4), JPTLPR(4), RSSQ(4), FRPRSQ(
     24), IDENO(4), JDERPR(4), LPLACE(4), P1(10,2,2)
      DIMENSION Q(4),Y(4),B(4),W(4),N1(4), P(4), G(4),DL(4),DK(4),UPD(4)
     1,UPB(4),UPT(4),T(8),     DB(4),DE(4),DY(4),ER(4)
      DIMENSION A(1,4,4),DD(4,4)
      K=SQRT(EF)
      X=RS(NATOM)**2
      NA=NATOM
      DO 11 I=1,4
      DD(I,NA)=DD(I,NA)/X
      A(1,I,NA)=A(1,I,NA)/X
      write(50,*) DD(I,NA)
      write(50,*) A(1,I,NA)
  11  continue
      Z=K*RS(NATOM)
      C=COS(Z)
      S=SIN(Z)
      U=RS(NATOM)/(2.*K)
      V=1./(2*K)
      F=-C
C
C          BESSEL FUNCTIONS OF THE FIRST KIND
      B(1)=S/Z
      B(2)=(S/Z**2)-(C/Z)
      B(3)=(((3/Z**3)-(1./Z))*S)-((3./Z**2)*C)
      B(4)=(((15./Z**4)-(6./Z**2))*S)-(((15./Z**3)-(1./Z))*C)
C
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE FIRST KIND WIT
C     RESPECT TO R.
      P(1)=K*((C/Z)-(S /Z**2))
      P(2)=K*((2*C/Z**2)-(2*S/Z**3)+(S/Z))
      P(3)=K*((9*C/Z**3)-(9*S/Z**4)+(4*S/Z**2)-(C/Z))
      P(4)=K*((60*C/Z**4)-(60*S/Z**5)-(7*C/Z**2)+(27*S/Z**3)-(S/Z))
C
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE FIRST KIND WIT
C     RESPECT TO E.
      DO  221 I=1,4
  221 DB(I)=(U/K)*P(I)
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE FIRST KIND WIT
C     RESPECT TO E AND R.
      DE(1)=V*(((1./(Z**2))-1.)*S-(1./Z)*C)
      DE(2)=V*(((4./(Z**3))-(2./Z))*S+(1.-(4./(Z**2)))*C)
      DE(3)=V*(((4./Z)-(27./(Z**3)))*C-((13./(Z**2))-(27./(Z**4))-1.)*S)
      DE(4)=V*(((240./(Z**5))-(114./(Z**3))+(7./Z))*S-((240./(Z**4))-(34
     1./(Z**2))+1.)*C)
C
C          BESSEL FUNCTIONS OF THE SECOND KIND.
      Y(1)=F/Z
      Y(2)=(F/Z**2)-(S/Z)
      Y(3)=(((3/Z**3)-(1./Z))*F)-((3./Z**2)*S)
      Y(4)=(((15./Z**4)-(6./Z**2))*F)-(((15./Z**3)-(1./Z))*S)
C
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE SECOND KIND WI
C     RESPECT TO R.
      G(1)=K*((S/Z)-(F /Z**2))
      G(2)=K*((2*S/Z**2)-(2*F/Z**3)+(F/Z))
      G(3)=K*((9*S/Z**3)-(9*F/Z**4)+(4*F/Z**2)-(S/Z))
      G(4)=K*((60*S/Z**4)-(60*F/Z**5)-(7*S/Z**2)+(27*F/Z**3)-(F/Z))
C
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE SECOND KIND WI
C     RESPECT TO E.
      DO 441 I=1,4
  441 DY(I)=(U/K)*G(I)
C
C          THE DERIVATIVES OF THE BESSEL FUNCTIONS OF THE SECOND KIND WI
C     RESPECT TO E AND R.
      ER(1)=V*(((1./(Z**2))-1.)*F-(1./Z)*S)
      ER(2)=V*(((4./(Z**3))-(2./Z))*F+(1.-(4./(Z**2)))*S)
      ER(3)=V*(((4./Z)-(27./(Z**3)))*S-((13./(Z**2))-(27./(Z**4))-1.)*F)
      ER(4)=V*(((240./(Z**5))-(114./(Z**3))+(7./Z))*F-((240./(Z**4))-(34
     1./(Z**2))+1.)*S)
C
C          THE PHASE SHIFT.
      DO 115 I=1,4
      DL(I)=(A(1,I,NA)*B(I)-P(I))
      DK(I)=(A(1,I,NA)*Y(I)-G(I))
      Q(I)=ATAN(DL(I)/DK(I))
  115 CONTINUE
C
C          THE DERIVATIVE OF THE PHASE SHIFT.
      DO 1111 I=1,4
      UPD(I)=(1./DK(I))
      UPB(I)=(DL(I)/(DK(I)**2))
      UPT(I)=(1./(1.+((DL(I)/DK(I))**2)))
      T(I)=A(1,I,NA)*DB(I)+B(I)*DD(I,NATOM)-DE(I)
      T(I+4)=A(1,I,NA)*DY(I)+Y(I)*DD(I,NA)-ER(I)
      W(I)=UPT(I)*(UPD(I)*T(I)-UPB(I)*T(I+4))
 1111 CONTINUE
      RETURN
      END
      SUBROUTINE RADIAL(XC)
      IMPLICIT REAL*8(A-H,O-Z)       
      REAL*8 N1       
      COMMON DELTA,    DDELTA,EF,N1,RW,RS ,RSQLND,DLD
      COMMON R,POT,ANTLH,SHIFT,RSSQ,FRPRSQ,P1,
     1ICHG,KMAX,NUMBER,KSTART,
     2IDENO,LPLACE,JDERPR,LMAX2,NATOM,NMAX
      DIMENSION R(800,4), POT(800,4), ICHG(10,4), KMAX(4), ANTLH(4), SHI
     1FT(4), RS(4), NUMBER(4), KSTART(10,4), JPTLPR(4), RSSQ(4), FRPRSQ(
     24), IDENO(4), JDERPR(4), LPLACE(4), P1(10,2,2)
      DIMENSION RSQLND(1,4,4),DLD(4,4),OVER(4,4)
      DIMENSION DELTA(4),DDELTA(4),N1(4)
    5 READ(7,195)   N3,B,LMAX2,C,D,N,N2
      NATOM=N
      write(51,*) NATOM
      KMAX(N)=N3
      ANTLH(N)=B
      SHIFT(N)=C
      RS(N)=D
      RSSQ(N)=D*D
      FRPRSQ(N)=12.5663706*D*D
      NUMBER(N)=N2
C     Change for 10 doublings 
      READ(7,201) (ICHG(K,N),K=1,9)
  201 FORMAT(9I4)
  245 FORMAT (2D20.13)
   30 DO 35 K=1,N3
   35 READ(7,245)    R(K,N),POT(K,N)
c  40 DO 45 K=1,LMAX2
c  45 READ(7,210)     KSTART(K,N),P1(K,1,N),P1(K,2,N)
      WRITE (8,215) N,N2
      WRITE (8,220) B,C,N,D
      WRITE (8,235)
      DO 55 K=1,10
      J=ICHG(K,N)
      IF (J) 55,55,50
   50 WRITE (8,240) R(J,N)
   55 CONTINUE
c     WRITE (8,230) (KSTART(K,N),P1(K,1,N),P1(K,2,N),K=1,LMAX2)
      WRITE (8,225)
      NTOP=AMIN0 (50,N3)
      DO 60 J=1,NTOP
   60 WRITE (8,250) (R(K,N),POT(K,N),K=J,N3,50)
      DO 65 K=1,N3
   65 POT(K,N)=POT(K,N)-C
      DO 70 K=1,N3
      L=N3-K+1 
      write(*,*) r(l,n),rs(n) 
      IF (R(L,N)-RS(N)) 80,75,70
   70 CONTINUE
   75 WRITE(8,260)
  260 FORMAT('0',' ERROR IN 70 LOOP')
   80 LPLACE(N)=L
      CALL INTG(XC)
C
  195 FORMAT (I4,F10.7,I2,2F10.6,I4,A6)
  200 FORMAT (10I4)
 205  FORMAT(F7.5,F15.5)
  210 FORMAT (I4,2E20.8)
  215 FORMAT(25HPOTENTIAL FOR ATOM NUMBER,I4///15H IDENTIFICATION,A6//)
  220 FORMAT (6HANTLH=,F12.8/6HSHIFT=,F12.8/4H RS(,I1,2H)=,F12.8//)
  225 FORMAT (1H1,8X,1HR,9X,9HPOTENTIAL,11X,1HR,9X,9HPOTENTIAL,11X,1HR,9
     1X,9HPOTENTIAL,11X,1HR,9X,9HPOTENTIAL//)
  230 FORMAT (//16H STARTING VALUES/38H   M       P1(M)               P1
     1(M+1)/(I4,2E20.8))
  235 FORMAT (19H CHANGE OF SCALE AT//)
  240 FORMAT (8H    R=  ,F10.4)
  250 FORMAT (4(F15.9,F15.7))
      RETURN
      END
      SUBROUTINE RELOGR(E,JRI,MXKAP,ALAG,ADX,DX1,BGX,RP2,HOC,MP,RLG,
     1 LWAVE)
      IMPLICIT REAL*8(A-H,O-Z)              
c     DOUBLE PRECISION DABS,DFLOAT
c     REAL*8 P,Q,PP,QP,UNP,WNP,UNP2,WNP2,H1,H2,H3,H4,H5,AH1,AH2,
c    $ AH3,AH4,AH5,TEST
      COMMON /FANG/ FREL(800,10),GREL(800,10)
      DIMENSION DX1(800)
      DIMENSION BGX(800),RLG(20),SXK(4),SXM(4),P(800),Q(800),PP(800)
     1 ,QP(800),RATFG(20,2)
      DIMENSION ALAG(6),RAD(800)
      C= 274.07204
      TEST= 1.D+5
      CIN=C*C
      CIN=1.0/CIN
 48   DO 95 K=1,MXKAP
      FLK= FLOAT(K*(K-1))
C   STARTING VALUES
      IF(C.GT.1.0E3) GO TO 9
C   RELATIVISTIC CASE
      U= 1.0 + FLK - HOC*HOC
      write (59,*) HOC,FLK,U,C
      U= (SQRT(U) - 1.0)/HOC
      write (59,*) HOC,FLK,U,C
    7 U= U*C
      write (59,*) HOC,FLK,U,C
      GO TO 8
C   NON - RELATIVISTIC STARTING VALUES
    9 IF(K.EQ.1) GO TO 10
      U= FLOAT(K-1)/ADX
      GO TO 8
   10 U= -BGX(1)*ADX/2.
  8   P(1)=1.D-20
      Q(1)= U*P(1)
      write (58,*) Q(1),U,BGX(1),ADX
      PP(1)= (1.0 + (E + BGX(1))*CIN)*Q(1) + P(1)/ADX
      QP(1)= -Q(1)/ADX + ((FLK/(ADX*ADX*(1.0 + (E+BGX(1))*CIN))) -
     $ (E + BGX(1))) *P(1)
   11 X= ADX
      N=1
  25  IK=0
      XC=X
      BGC=BGX(N)
      WC=Q(N)
      UC=P(N)
  20  IK=IK+1
   12 SXM(IK)= ADX*( -WC/XC + ((FLK/(XC*XC*(1.0 + CIN*(E+BGC)))) -
     $ (E+BGC))*UC)
      SXK(IK)= ADX*((1.0 + (E+BGC)*CIN)*WC + UC/XC)
      write (57,*) SXM(IK),SXK(IK),ADX,WC,XC,FLK,CIN,E,BGC,UC
  15  GO TO (16,17,18,19),IK
   16 XC= XC+ .5*ADX
      UC=UC+.5*SXK(1)
      WC=WC+.5*SXM(1)
      BGC=.5*(BGC+BGX(N+1))
      GO TO 20
  17  UC=UC+.5*(SXK(2)-SXK(1))
      WC=WC+.5*(SXM(2)-SXM(1))
      GO TO 20
   18 XC= XC+ .5*ADX
      UC=UC+SXK(3)-.5*SXK(2)
      WC=WC+SXM(3)-.5*SXM(2)
      BGC=BGX(N+1)
      GO TO 20
  19  Q(N+1)=Q(N)+(SXM(1)+2.*SXM(2)+2.*SXM(3)+SXM(4))/6.
      P(N+1)=P(N)+(SXK(1)+2.*SXK(2)+2.*SXK(3)+SXK(4))/6.
      write (56,*) Q(N+1),P(N+1),Q(N),P(N)
      write (56,*) SXM(1),SXM(2),SXM(3),SXM(4)
      write (56,*) SXK(1),SXK(2),SXK(3),SXK(4)
      QP(N+1)= -Q(N+1)/XC + ((FLK/(XC*XC*(1.0 + (E+BGC)*CIN))) -
     $ (E+BGC))*P(N+1)
      PP(N+1)= (1.0 + (E+BGC)*CIN)*Q(N+1) + P(N+1)/XC
   24 X= X+ADX
      N=N+1
      MESH= N
      IF (N-6)25,26,26
   26 DX= DX1(N)
       X= X + DX
      RAD(N + 1)= X
      H1= 3.3D0
      H2= -4.2 D0
      H3= 7.8 D0
      H4= -4.2 D0
      H5= 3.3 D0
      AH1= .311111111111 D0
      AH2= 1.42222222222     D0
      AH3= .533333333333 D0
      AH4= 1.42222222222 D0
      AH5= .311111111111 D0
C     IF(N.GT.175) GO TO 650        
C     Change for 10 doublings
      IF(N.GT.303) GO TO 650                                            
      IF(MESH.NE.32) GO TO 600
      H1= 6.32430555556 D0
      H2= -13.4847222222 D0
      H3= 16.2166666667 D0
      H4= -8.52638888889 D0
      H5= 2.97013888889 D0
      AH1= .295138888889 D0
      AH2= 1.73611111111 D0
      AH3= -.868055555556D0
      AH4= 1.30208333333 D0
      AH5= .034722222222 D0
      GO TO 650
  600 CONTINUE
      IF(MESH.NE.33) GO TO 601
      H1= 3.66222222222 D0
      H2= -11.2888888889 D0
      H3= 19.9111111111 D0
      H4= -12.2666666667 D0
      H5= 3.98222222222 D0
      AH1= .315 D0
      AH2= 1.425 D0
      AH3= .225 D0
      AH4= .96 D0
      AH5= .075 D0
      GO TO 650
  601 CONTINUE
      IF(MESH.NE.34) GO TO 602
      H1= 3.1865625 D0
      H2= -3.9796875     D0
      H3= 10.9546875 D0
      H4= -9.36 D0
      H5= 3.6984375 D0
      AH1= .329340277778 D0
      AH2= 1.32197916667 D0
      AH3= .774131944444 D0
      AH4= 1.01232638889 D0
      AH5=.0622222222222 D0
      GO TO 650
  602 CONTINUE
      IF(MESH.NE.35) GO TO 603
      H1= 3.04563492063 D0
      H2= -2.95833333333     D0
      H3= 5.48611111111 D0
      H4= -2.98611111111 D0
      H5= 2.41269841270 D0
      AH1= .311111111111 D0
      AH2= 1.42222222222 D0
      AH3= .533333333333 D0
      AH4= 1.42222222222 D0
      AH5= .311111111111 D0
      GO TO 650
  603 CONTINUE
      IF(MESH.NE.36) GO TO 650
      H1= 3.00590277778 D0
      H2= -2.71371527778 D0
      H3= 4.77239583333 D0
      H4= -1.05225694444 D0
      H5=1.48767361111 D0
      AH1= .311111111111 D0
      AH2= 1.42222222222 D0
      AH3= .533333333333 D0
      AH4= 1.42222222222 D0
      AH5= .311111111111 D0
      MESH= 4
  650 CONTINUE
   27 UNP= P(N-5) + DX*(H1*PP(N)+ H2*PP(N-1)+ H3*PP(N-2)+ H4*PP(N-3)
     1 +H5*PP(N-4))
      WNP= Q(N-5) + DX*(H1*QP(N)+ H2*QP(N-1)+ H3*QP(N-2)+ H4*QP(N-3)
     1 + H5*QP(N-4))
      NIT=0
   33 QP(N+1)= -WNP/X + ((FLK/(X*X*(1.0 + CIN*(E+BGX(N+1))))) -
     $ (E+ BGX(N+1)))*UNP
      PP(N+1)= (1.0 + CIN*(E + BGX(N+1)))*WNP + UNP/X
      UNP2= P(N-3)+ DX*(AH1*PP(N+1)+ AH2*PP(N)+ AH3*PP(N-1)+ AH4*PP(N-2)
     1 + AH5*PP(N-3))
      WNP2= Q(N-3)+ DX*(AH1*QP(N+1)+ AH2*QP(N)+ AH3*QP(N-1)+ AH4*QP(N-2)
     1 + AH5*QP(N-3))
      IF(DABS(TEST*(UNP2-UNP))-DABS(UNP2))30,30,31
  30  IF(DABS(TEST*(WNP2-WNP))-DABS(WNP2))32,32,31
   31 IF(NIT-20) 81,32,81
  81  NIT=NIT+1
      WNP=WNP2
      UNP=UNP2
      GO TO 33
   32 Q(N+1)= WNP2
      P(N+1)= UNP2
      N=N+1
      MESH= MESH + 1
      IF (N-JRI)26,35,35
   35 IF(LWAVE) 34,34,36
   36 DO 37  MWAVE=1,JRI
      FREL(MWAVE,K)= Q(MWAVE)/C
      GREL(MWAVE,K)= P(MWAVE)
      write (55,*) FREL(MWAVE,K),GREL(MWAVE,K),Q(MWAVE),P(MWAVE)
   37 CONTINUE
   34 CONTINUE
      SOP1= .0 D0
      SOP2= .0 D0
      DO 655  KK=1,6
      L= KK + MP
      P(L)= P(L)/RAD(L)
      Q(L)= Q(L)/RAD(L)
      write (55,*) P(L),Q(L),RAD(L)
      SOP1= SOP1+ ALAG(KK)* P(L)
      SOP2= SOP2 + ALAG(KK)*((1.0 + CIN*(E+ BGX(L)))*Q(L))
      write (55,*) SOP1,SOP2,RP2
  655 CONTINUE
      RLG(K)= RP2*SOP2/SOP1
      write (55,*) RLG(K)
  95  CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIMCHG (F,RR,N,LPLACE,RS)
      IMPLICIT REAL*8(A-H,O-Z)              
      DIMENSION LPLACE(4), RS(4)
      DIMENSION F(800), RR(800)
      X=RR(1)*(4.*F(1)+F(2))
      IMAX=LPLACE(N)+1
      DO 5 I=4,IMAX,2
      H=RR(I)-RR(I-1)
    5 X=X+H*(F(I)+F(I-2)+4.*F(I-1))
      I=IMAX- MOD (IMAX,2)
      EPS=RS(N)-RR(I)
      SIMCHG=X/3.+EPS*F(I)+EPS*EPS*(F(IMAX)-F(IMAX-1))/(2.*H)
      RETURN
      END
