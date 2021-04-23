       PROGRAM INTPOT
c      this program interpolates lapw potentials to the apw mesh.
c     DX=0.015
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(1001),RV(1001),REE(1001),VNEW(1001),B(1001)
      DIMENSION Y(1001),V(1001)
      DIMENSION TITLE(10)
      REAL*8 INTERP
      OPEN(5,FILE='input.dat',BLANK='ZERO',STATUS='UNKNOWN')
      OPEN(7,FILE='apwpot.dat',BLANK='ZERO',STATUS='UNKNOWN')
      OPEN(10,FILE='lapwpot.dat',BLANK='ZERO',STATUS='UNKNOWN')
c     NHS=368
c     VT=0.0
c      z is not used in this program
c     Z=56.0
      PI=3.141593
      LMAX=10
      READ (5,1493)TITLE,IATOMS
      write(6,1494)TITLE,iatoms
1493  FORMAT (10A4,I5)
1494  FORMAT (10A4,i5)
1495  FORMAT (10A4)
      read(10,1495) skip1
      read(10,1495) skip1
      read(10,1495) skip1
      read(10,1495) skip1
      read(10,1495) skip1
      DO 5000 IATOM=1,IATOMS
      DO 5001 IN=1,1001
      R(IN)=0.0
      RV(IN)=0.0
      REE(IN)=0.0
      VNEW(IN)=0.0
      B(IN)=0.0
      Y(IN)=0.0
 5001 CONTINUE 
      READ (5,160)TITLES
 160  FORMAT (A2)
      READ (5,140)NHS,VT
      write (6,140)NHS,VT
 140  FORMAT (I5,F10.5)
c     WRITE (7,159)
      READ (5,150) N,L,IDOUB,RS,ISPACE
      write(6,150) N,L,IDOUB,RS,ISPACE
150   FORMAT (3I5, F10.6,I8)
      NN=N+1
      NNN=NHS+1
 164  FORMAT (F8.6)
c     DX=0.015
      DX=0.03
C     RNOT=0.15073307E-03
C     WRITE(6,9015)1,RMUF,RNOT
c9015 FORMAT(' MUFFIN TIN RADIUS FOR ATOM #',I1,' IS',
c    *        E16.9,/' RNOT IS',1X,E16.9)
 9016 FORMAT(1X,'MUFFIN TIN ZERO IS=',F10.6)
C     IF (ISRUC.EQ.2) GO TO 130
C     IF (ISRUC.EQ.1) GO TO 120
 120  CONTINUE
      if (iatom.eq.1) read(10,*) skip1,skip2,skip3,skip4,skip5
      if (iatom.eq.1) read(10,*) skip1,skip2,skip3
      if (iatom.eq.2) read(10,*) skip1,skip2,skip3,skip4,skip5
      if (iatom.eq.2) read(10,*) skip1,skip2,skip3
      if (iatom.eq.2) read(10,*) skip1,skip2,skip3
      if (iatom.eq.2) read(10,*) skip1,skip2,skip3
      DO 9004 J=1,NHS
9004  R(J)=RS*EXP(-DX*(NHS-J))      
C     READ 1,R(I), V(I)
C   1 FORMAT (F9.5,F15.6)
      READ (10,1)(RV(I),I=1,NHS)
C     DO 10 I=1,NHS
C     RV(I)=-RV(I)/SQRT(4.*PI)
c     WRITE(6,123) R(j),RV(j)
    1 FORMAT(4F20.0)
c   1 FORMAT (3X,D21.14,4X,D21.14,4X,D21.14,4X,D21.14)
  123 FORMAT (2E15.6)
C   1 FORMAT (F9.5,E11.5)
   10 CONTINUE
      DO 12 II=1,NHS
      WRITE(6,123) R(II),RV(II)
 12    CONTINUE
C     R(1)=0.0
C     RV(1)=2.0*Z
      write(6,250) RS
 250  FORMAT (F15.6,'rs')
 260  H=RS*1.0E5/ISPACE
      H=AINT(H+0.5)/1.0E6
c     write(6, 999) H
 999  FORMAT (1X,F15.6,'h')
      NP=NN+1
       H=0.000025
c      H=0.00002
      CALL INDEXX (REE,NP,L,IDOUB,H)
      HH=REE(NN)-REE(N)
      CHECK=RS+6.*HH
      IF (REE(NN).LT.CHECK) ISPACE=ISPACE-10
      IF (REE(NN).LT.CHECK)  GO TO 260
      J=2
      N2=N-1
      DO 45 K=2,N2
   25 IF (REE(K)-R(J))   40,35,30
   30 J=J+1
      IF (J.LT.NNN) GO TO 25
      GO TO 300
 35   B(K)=RV(J)
      VNEW(K)=B(K)/REE(K)
      GO TO 45
   40   B(K)      =INTERP(REE(K),R,RV(1),J,3)
      VNEW(K)=B(K)/REE(K)
      WRITE( 6,200) REE(K ),VNEW(K)
   45 CONTINUE
 300  NRSS=K-1
      write(6,301) NRSS
 301  FORMAT (1X,I4,'nrss')
      Y(NRSS-1)=VNEW(NRSS-1)  *REE(NRSS-1)
      Y(NRSS)=VNEW(NRSS)  *REE(NRSS)
      DO 11 I=NRSS,NN
      Y(I+1)=2.*Y(I)-Y(I-1)
      VNEW(I+1)=  Y(I+1)/REE(I+1)
  11  CONTINUE
      ICYCLE=1
      IATO=IATOM
C     WRITE (7,139)
 139  FORMAT (4HPOTE)
C     WRITE (6,155) N,     REE(2),LMAX,VT,RS,IATO,TITLES,   ICYCLE
      WRITE (7,155) N,     REE(2),LMAX,VT,RS,IATO,TITLES,   ICYCLE
 155  FORMAT (I4,F10.7,I2,2F10.6,I4,A2,I2)
C     WRITE (6,165)
      WRITE (7,165)
c165  FORMAT (20H  32  64  96 128 160)
165   FORMAT (36H  32  64  96 128 160 192 224 256 288)
      
  201 FORMAT (2D20.13)
C165  FORMAT (16H  32  64  96 128)
C     WRITE( 6,200)(REE(JJ),VNEW(JJ),JJ=2,NN)
      WRITE (7,201)(REE(JJ),VNEW(JJ),JJ=2,NN)
c201  FORMAT (F8.5,F15.8)
 200  FORMAT (1X,F7.5,F15.6)
  135 FORMAT (4H   1,2E20.8)
    4 FORMAT(1H ,F7.5,F15.6)
    5 FORMAT(F7.5,F15.6)
 5000 CONTINUE
      STOP
      END
      SUBROUTINE INDEXX (REE,N1,L,IDOUB,H)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION REE(N1)
      I=1
      M=IDOUB
      DELTAX=H
      REE(1)=0.0
      DO 10 J=1,L
      IF (L.EQ.J) M=N1-IDOUB*(L-1)-1
      DO 5 K=1,M
      I=I+1
    5 REE(I)=REE(I-1)+DELTAX
      DELTAX=DELTAX+DELTAX
   10 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION INTERP(XX,X,F,J1,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(1),F(1)
      FX=0.
      ISTART=J1-N/2
      J2=ISTART+N-1
      IF(ISTART.GE.1) GO TO 60
      J2=J2+1-ISTART
      ISTART=1
   60 CONTINUE
      DO 10 J=ISTART,J2
      P=F(J)
      DO 5 I=ISTART,J2
      IF (I.EQ.J) GO TO 5
      P=P*(XX-X(I))/(X(J)-X(I))
    5 CONTINUE
   10 FX=FX+P
      INTERP=FX
      RETURN
      END
