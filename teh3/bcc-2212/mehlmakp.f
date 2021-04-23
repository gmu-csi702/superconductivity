      PROGRAM EXBIRPLT
C
C     Uses the extended Birch energy-volume parameterization,
C
C     E(v)=Sum[a(n) v^(-2n/3),{n,0,N}]
C
C     to calculate the equation of state,
C
C     P(v)=v^(-5/3) Sum[(n+1) a(n+1) v^(-2n/3),{n,0,N-1}]
C
C     Output is either E(V), P(V), V(P) or G(P), and pressure is expressed
C      in Mbar and volumes in cubic Angstroms.  Volume dependence may be
C      changed to lattice constant dependence for several lattices.
C
C     Added the ability to calculate B(V) or B(P).  -mjm  10-Dec-1986
C
C     Bonus!  We can print out distances in au or Angstroms
C      & energies in eV or Ry!                      -mjm  20-Sep-1989
C
C     The polynomial coefficients a(n) and the equilibrium volume
C       VO are created by the fitting program EXBIRFT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXO=25)
      COMMON A(0:MAXO),EARRAY(0:1)
      CHARACTER*70 POLYFL,EOSFL
      CHARACTER*1 BLANK
      CHARACTER*30 VSTRING(4)
C
C     Data below is 1 Ry in ev, 1 au^3 in Angstrom^3, and
C       1 Ry/au^3 in Mbar, using the 1986 determinations of these
C       constants
C
      DATA RY/1.36056981D1/,AU3/1.481847435D-1/,
     1     AU/5.29177249D-1/,
     2     RYAU3/1.47105164D2/,TWO3/6.666666666666666667D-1/
      DATA ZERO/0D0/,ONE/1D0/
      DATA BLANK/' '/
      DATA VSTRING/'volume in Angstrom^3:         ',
     1             'lattice constant in Angstroms:',
     2             'volume in au^3:               ',
     3             'lattice constant in au:       '/
C
      OPEN (5,FILE='plot5.dat',BLANK='ZERO',STATUS='UNKNOWN')
      READ(5,27) POLYFL
27    FORMAT(A70)
      IF(POLYFL.EQ.BLANK) STOP
      OPEN(UNIT=11,FILE=POLYFL,STATUS='OLD',ERR=1400)
C
C     Get the parameter information out of this file:
C
      READ(11,*) N
      READ(11,*) (A(I),I=0,N)
      CLOSE(11)
      WRITE(6,135) N,(A(I),I=0,N)
135   FORMAT(/1X,I5,1P/(1X,5E15.7))
C   IVOPT=1 for sc, IVOPT=2 for fcc, IVOPT=3 for bcc
      READ(5,136) IOPT,IENG,IVOPT,IUDIS
136   FORMAT (I5,I5,I5,I5)
C
C     Default is P(V)
C
      IF(IOPT.LT.1.OR.IOPT.GT.6) IOPT=2
C     If IOPT=1, how should the energy be printed?
      IF(IOPT.EQ.1) THEN
       IF(IENG.NE.2) THEN
        ESCALE=RY
       ELSE
        ESCALE=ONE
       END IF
      ELSE
C      Leave it in Ry:
       ESCALE=ONE
      END IF
C
C     If a volume dependence is asked for, allow it to be replaced
C      by a lattice constant dependence in certain situations:
      IF((IOPT.EQ.1).OR.(IOPT.EQ.2).OR.(IOPT.EQ.5)) THEN
C      Angstroms or au:
       IF(IUDIS.NE.2) IUDIS=1
C      Volume and length scale factors:
       IF(IUDIS.EQ.1) THEN
        VOLSCL=ONE/AU3
        DISSCL=ONE/AU
       ELSE
        VOLSCL=ONE
        DISSCL=ONE
       END IF
C      Scale factor for volume determination (V=SCALE*A*A*A):
       IF(IUDIS.EQ.1) THEN
        IVOPT2=2
       ELSE
        IVOPT2=4
       END IF       
       IF(IVOPT.EQ.1) THEN
        SCALE=1D0
       ELSE IF(IVOPT.EQ.2) THEN
        SCALE=2.5D-1
       ELSE IF(IVOPT.EQ.3) THEN
        SCALE=5D-1
       ELSE
C       Pure volume is the default:
        IVOPT=4
        IVOPT2=IVOPT2-1
       END IF
      END IF
C
      READ(5,27) EOSFL
      OPEN(UNIT=21,FILE=EOSFL,STATUS='NEW')
C
      READ(5,195) XMIN,XMAX,DELX
C     WRITE(6,195) XMIN,XMAX,DELX
195   FORMAT (F8.3,F9.3,F7.3)
      DO 200 XI=XMIN,XMAX,DELX
C      Convert to volume if necessary:
196   FORMAT (F5.5)
      IF(IVOPT.LT.4) THEN 
          V=SCALE*XI*XI*XI
      ELSE
        V=XI
       END IF
C
C      Convert volume to au^3 (if necessary) and find polynomial parameter
       VAU=V*VOLSCL
       X=VAU**(-TWO3)
C
C      Call the polynomial to evaluate the energy and pressure:
C
       CALL PLYEVL(N,ZERO,A,X,2,EARRAY)
C
C      Energy has been calculated in Ry.  If that's all we want, print
C       it out (in eV)
       IF(IOPT.EQ.1) THEN
C       Note that XI has the appropriate units:
        FIRST=XI
        SECOND=ESCALE*EARRAY(0)
       ELSE
C       We'll need the pressure:
        P=TWO3*X*EARRAY(1)/VAU
        IF(IOPT.EQ.2) THEN
C        P(V) in Mbar
         FIRST=XI
         SECOND=RYAU3*P
        ELSE IF(IOPT.EQ.3) THEN
C        V(P)
         FIRST=RYAU3*P
         SECOND=V
        ELSE IF(IOPT.EQ.4) THEN
C        G(P), in eV as a function of P in Mbar
         FIRST=RYAU3*P
         SECOND=RY*(EARRAY(0)+P*VAU)
        ELSE
C        We'll need the bulk modulus:
         BULKM=RYAU3*(X**(2.5D0))*(4D0*X*EARRAY(2)+
     1         1D1*EARRAY(1))/9D0
         IF(IOPT.EQ.5) THEN
          FIRST=XI
          SECOND=BULKM
         ELSE IF (IOPT.EQ.6) THEN
          FIRST=RYAU3*P
          SECOND=BULKM
         END IF
        END IF
       END IF
C
C      Output information into file:
C
200    WRITE(21,205) FIRST,SECOND
205    FORMAT(1X,2F18.8)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       STOP
1400  WRITE(*,1405)
1405  FORMAT ('ERROR IN INPUT')
      END
C________________________________________________________________________
C
      SUBROUTINE  P L Y E V L (M,X0,A,X,N,P)
C     Evaluates the polynomial given by
C
C       p(x) = Sum [ a(i) (x-x0)^i ,{i,0,m}]
C
C       and its first N derivatives
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Note the dummy indexing on A and P:
      DIMENSION A(0:1),P(0:1)
      DATA ZERO/0D0/
      Y=X-X0
C     Zeroth order term (in Y)
      IPROD=1
      DO 100 J=0,N
        P(J)=IPROD*A(M)
100     IPROD=IPROD*(M-J)
      DO 200 I=M-1,0,-1
        IPROD=1
        DO 200 J=0,N
          IF(IPROD.GT.ZERO) P(J)=P(J)*Y+IPROD*A(I)
200       IPROD=IPROD*(I-J)
      RETURN
      END
