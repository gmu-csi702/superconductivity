      PROGRAM EXBIRFT
C     Fits energy versus volume data to the form
C
C     E(V) = Sum [ a(n) v^(-2n/3) , {n,0,N}]
C
C     and then calculates E(Vo), Vo, Ko, and Ko'
C
C     This version uses the Least Squares fitting routine of
C      the Numerical Recipes book.
C
      IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)
      PARAMETER (MAXN=1000,MAXM=10,MAXM1=MAXM+1)
      PARAMETER (MAXCOL=20)
      LOGICAL CONA,DISB,ENGR,MINFOUND
      COMMON V(MAXN),E(MAXN),A(0:MAXM),P(0:MAXM,MAXN),
     1   EARRAY(0:3),SIG(MAXN),CVM(0:MAXM,0:MAXM),ASIG(0:MAXM),
     2   COL(MAXCOL),LISTA(0:MAXM)
C     CONV can't be in COMMON because Professional Fortran doesn't
C      like DATA initialization of COMMON:
      DIMENSION CONV(3)
      CHARACTER*50 DATAFL,OUTFL
      CHARACTER*75 TITLE
      CHARACTER*1 BLANK,ANS,YESU,YESL
C     Conversion factors:  V=CONV*(A**3)
      DATA CONV/2.5D-1,1D0,5D-1/
      DATA BLANK/' '/,YESU/'Y'/,YESL/'y'/
C     New (1986) Values for 1 Bohr in  and 1 Rydberg in eV:
      DATA AUOAN/5.29177249D-1/,RYOEV/1.36056981D1/
      DATA THIRD/3.33333333333333333D-1/,TWO3/6.666666666666666667D-1/
      DATA ONE/1D0/
      OPEN (5,FILE='fit5.dat',BLANK='ZERO',STATUS='UNKNOWN')

C
C     We'll assume all fitting parameters need to be fitted.  Note that
C      LISTA requires 1 + our index:
      DO 10 I=0,MAXM
10     LISTA(I)=I+1
C
      AU3=AUOAN*AUOAN*AUOAN
      READ(5,115,ERR=1000,END=1000) DATAFL
115   FORMAT(A50)
      IF(DATAFL.EQ.BLANK) STOP
      OPEN(UNIT=22,FILE=DATAFL,STATUS='OLD',ERR=998)
C   LAT=1 for fcc LAT=2 for sc LAT=3 for bcc 
      READ(5,116) IWHICH,LAT,IAUNIT,IEUNIT
116   FORMAT (I5,I5,I5,I5)
      CONA=IWHICH.EQ.2
      IF(.NOT.CONA) IWHICH=1
C! Default is volume
      DISB=IAUNIT.EQ.2
C     Default is Bohrs:
      IF(.NOT.DISB) IAUNIT=1
      IF(DISB) THEN
       AUNIT=1D0
      ELSE
       AUNIT=AUOAN
      END IF
      VUNIT=AUNIT*AUNIT*AUNIT
      ENGR=IEUNIT.EQ.2
C     Default is Rydbergs:
      IF(.NOT.ENGR) IEUNIT=1
      IF(ENGR) THEN
       EUNIT=1D0
      ELSE
       EUNIT=RYOEV
      END IF
C
C     To find the minimum energy, we'll need an estimated starting value.
C      Since we've got the energies on file, we might as well use the
C      volume of the lowest energy as our estimated Vo.
      EMIN=1D10
      VMIN=0D0
C     Read data in the form X,E (x=volume or lattice constant, E=energy)
C     where Angstroms and eV or atomic units (Bohrs and Rydbergs) are
C     used depending upont the setting of IUNIT 
C   The parameter "M" is the order of the polynomial.
 
      READ (5,185) LSKIP,ICOL,ICX,ICE,M
185   FORMAT (I5,I5,I5,I5,I5)
      DO 190 LS=1,LSKIP
       READ(22,191) TITLE
190    WRITE(*,191) TITLE
191    FORMAT(A75)
      N=0
200   READ(22,*,ERR=200,END=300) (COL(I),I=1,ICOL)
      N=N+1
      X=COL(ICX)
      EE=COL(ICE)
      IF(CONA) THEN
       VA=CONV(LAT)*X*X*X
      ELSE
       VA=X
      END IF
C     Volume in au**3:
      V(N)=VA/VUNIT
C     Energy in Rydbergs:
      E(N)=EE/EUNIT
C     For now assume that all data are created equal:
      SIG(N)=ONE
      PRINT 255, N,VA,EE,V(N),E(N)
255   FORMAT(1X,I5,4F15.5)
      IF(E(N).LT.EMIN) THEN
       EMIN=E(N)
       VMIN=V(N)
      END IF
      GO TO 200
300   CLOSE(22)
C     XSTART is the starting value for the minimum search:
      XSTART=VMIN**(-TWO3)
      PRINT 305, EMIN,VMIN,XSTART
305   FORMAT(/' Minimum energy in file = ',F15.5,' at V = ',
     1        2F15.7/)
C
C     Set up the fit.  Note that the functions to be fitted are
C      v^(-2M/3),n=0,1,2,...MAXM
      M1=M+1
      DO 400 I=1,N
C      Establish the basis functions:
       P(0,I)=1D0
       X=V(I)**(-TWO3)
       DO 400 J=1,M
400     P(J,I)=X**J
      CALL LFIT(E,SIG,N,A,M1,LISTA,M1,CVM,MAXM1,CHISQ,P)
C     Figure the rms "error" in each term:
C     CHI is the numerical recipies definition, unless N=2, in which
C      case we won't get a good fit anyway
      IF(N.GT.2) THEN
       CHI=SQRT(CHISQ/(N-2))
      ELSE
       WRITE(*,414)
414    FORMAT(/' Warning:  N<3, so delta A''s are not well defined')
       CHI=SQRT(CHISQ)
      END IF
      DO 410 I=0,M
410    ASIG(I)=SQRT(CVM(I,I))*CHI
      PRINT 415, (I,A(I),ASIG(I),I=0,M)
415   FORMAT(/' Fitting coefficients:'/(1X,I5,1P2E16.8))
C
C     Now for the error checking:
C
      ERMS=0D0
      EMAX=0D0
      DO 600 I=1,N
       XO=V(I)**(-TWO3)
       CALL PLYEVL(M,0D0,A,XO,0,EARRAY)
       ECK=EARRAY(0)
       ERR=ECK-E(I)
       IF(ABS(ERR).GT.ABS(EMAX)) EMAX=ERR
       ERMS=ERMS+ERR*ERR
600    PRINT 605, I,V(I),E(I),ECK,ERR
605    FORMAT(1X,I5,F12.5,3F15.5)
      ERMS=SQRT(ERMS/N)
C     Convert the errors to eV
      ERMSEV=ERMS*RYOEV
      EMAXEV=EMAX*RYOEV
C     Now we must find the equilibrium volume, VO, if we can.  We'll
C      use Newton's method.  Note that we can write the energy
C      as E(v)=f(v^(-2/3)).  Then every extremum of E is also an
C      extremum of f.
C     Pick a trial value for XO.  Eventually, VO=XO**(-3/2).
C      The logical starting value is VMIN**(-2/3), which we
C      calculated above.  After the first minimization, we should
C      have a good estimate of VO, so we'll use that.
      XO=XSTART
C     Since we can calculate exact second derivatives of the
C      fit for M>1, we can use Newton's method to find
C      the root.  Try no more than 100 times
      DO 700 IT=1,100
       CALL PLYEVL(M,0D0,A,XO,2,EARRAY)
C      EARRAY contains f and its derivatives
       DIFF=EARRAY(1)/EARRAY(2)
       PRINT 695, XO,EARRAY(1),EARRAY(2),DIFF
695    FORMAT(1X,1P4E16.7)
       XO=XO-DIFF
       IF(ABS(DIFF).LT.1D-5) GO TO 800
700    CONTINUE
C     No roots found.  Ask if fit should be printed out anyway:
      PRINT 705
705   FORMAT(/' No minimum found after 100 iterations.'
     1       /' Print out fitting parameters?  ')
      READ(5,875) ANS
      MINFOUND=.FALSE.
      GO TO 900
C
C     Use this new value of XO as the new starting value
800   XSTART=XO
      MINFOUND=.TRUE.
C     XOI is the 2/3rd root of the volume:
      XOI=1D0/XO
      VO=XOI**1.5D0
C     Now for the other equilibrium constants:
      CALL PLYEVL(M,0D0,A,XO,3,EARRAY)
      PRINT 711, XO,(EARRAY(I),I=0,3)
711   FORMAT(1X,1P5E15.6)
      EO=EARRAY(0)
C     To compare with "universal" equations of state, find the
C      predicted energy at infinity:
      EINF=A(0)
      ECOH=EINF-EO
C     Check that we've really found an extremum:
      PO=-TWO3*EARRAY(1)*(XO**2.5D0)
      PRINT 715, PO
715   FORMAT(' "Equilibrium" pressure = ',1PE15.6)
      KO=(1D1*EARRAY(1)+4D0*XO*EARRAY(2))*XO/(9D0*VO)
      KOP=(2.5D1*EARRAY(1)+4D0*XO*(6D0*EARRAY(2)+XO*EARRAY(3)))/
     1    (3D0*(5D0*EARRAY(1)+2D0*XO*EARRAY(2)))
      ALAT=(VO/CONV(LAT))**THIRD
C
C     Now use the ASIG to estimate the error terms.  We'll ignore
C      the correlation between the errors for now, and assume that
C      all terms are additive:
C
C     We'll only calculate DELVO, DELEO, and DELKO.  All involve
C      polynomials in XO
      DELEO=ZERO
      DELVO=ZERO
      DELKO=ZERO
      DO 850 I=N,0,-1
       ABSAG=ABS(ASIG(I))
       DELEO=DELEO*XO+ABSAG
       TX=TWO3*I
       DELVO=DELVO*XO+ABSAG*TX
850    DELKO=DELKO*XO+ABSAG*TX*(TX+ONE)
C     DELVO and DELKO need corrections:
      DELVO=DELVO/KO
      DELKO=DELKO/VO+ABS(KO*KOP*DELVO/VO)
      DELAO=THIRD*DELVO*ALAT/VO
      PRINT 835, M,VO,ALAT,EO,ECOH,KO,KOP,
     1           ERMS,EMAX
835   FORMAT(///' Equilibrium parameters for the Birch-Murnaghan',
     1 ' equation of order ',I2,':'/'    Vo = ',F15.5,
     2 ' Bohr**3'/5X,'a = ',F15.5,' Bohrs',/'    Eo = ',
     3 F15.5,' Rydbergs'/'    Ec = ',F15.5,' Rydbergs'
     4 /'    Ko = ',F15.5,' Rydbergs/Bohr**3'/'    Ko''= ',
     5 F15.5/'    RMS error in energy fit = ',F13.5,' Rydbergs'
     6 /'    Maximum error in energy fit = ',F9.5,' Rydbergs')
C     Now convert to standard units:
C     Energies in eV:
      EO=EO*RYOEV
c     DELEO=DELEO*RYOEV
c     ECOH=ECOH*RYOEV
C     Pressure in Mbar (1986 conversion factor)
      KO=1.47105164D2*KO
c     DELKO=1.47105164D2*DELKO
C     Volume in ^3
C     VO=VO*AU3
C     DELVO=DELVO*AU3
C     Lattice constant in 
c     ALAT=ALAT*AUOAN
c     DELAO=DELAO*AUOAN
      PRINT 855, M,VO,ALAT,EO,ECOH,KO,KOP,
     1           ERMSEV,EMAXEV
855   FORMAT(//' Equilibrium parameters for the Birch-Murnaghan',
     1 ' equation of order ',I2,':'/'    Vo = ',F15.5,
     2 ' Bohr**3'/5X,'a = ',F15.5,' BOHR',
     3/'    Eo = ',F15.5,'eV'/'    Ec = ',F15.5,' eV'
     4/'    Ko = ',F15.5,' Mbar'/'    Ko''= ',F15.5
     5/'    RMS error in energy fit = ',F13.5,' eV'
     6/'    Maximum error in energy fit = ',F9.5,' eV'//)
      READ(5,875) ANS
875   FORMAT(A1)
900   IF((ANS.EQ.YESL).OR.(ANS.EQ.YESU)) THEN
       READ(5,115) OUTFL
       OPEN(UNIT=23,FILE=OUTFL,STATUS='NEW')
       WRITE(23,895) M,(A(I),I=0,M)
895    FORMAT(1X,I5/(1P4E20.12))
C      Write out the estimated RMS variances:
       WRITE(23,8951) (ASIG(I),I=0,M)
8951   FORMAT(' RMS variances:'/(1P4E20.12))
C      Add the equilibrium functions as a reference:
       IF(MINFOUND) WRITE(23,897) EO,ECOH,VO,ALAT,
     1     KO,KOP
897    FORMAT(' Eo = ',F15.5,' eV'/' Ec = ',F15.5/
     1  ' Vo = ',F15.5,' Bohr^3'/' ao = ',F15.5,
     2  ' Bohr'/' Ko = ',F15.5,' MBar'
     3  /' Ko''= ',F15.5)
C      Also print the errors in the calculation
       WRITE(23,925) ERMSEV,EMAXEV
925    FORMAT(' RMS error = ',1PE11.3,' eV'/' Max error = ',
     1   E11.3,' eV')
       CLOSE(23)
      ENDIF
      STOP
 998  WRITE (6,999)
 999  FORMAT ('ERROR IN INPUT FILE,PROGRAM STOPPED')
1000  STOP
      END
C_______________________________________________________________________
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
      Y=X-X0
C     Zeroth order term (in Y)
      IPROD=1
      DO 100 J=0,N
        P(J)=IPROD*A(M)
100     IPROD=IPROD*(M-J)
      DO 200 I=M-1,0,-1
        IPROD=1
        DO 200 J=0,N
          IF(IPROD.GT.0D0) P(J)=P(J)*Y+IPROD*A(I)
200       IPROD=IPROD*(I-J)
      RETURN
      END
C_______________________________________________________________________
C
C     Numerical Recipes routines.  See pages 513 ff of the Fortran edition
C
C_______________________________________________________________________
C
      SUBROUTINE LFIT(Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=1000,MAXM=10,MAXM1=MAXM+1)
      PARAMETER (MMAX=MAXM1)
      DIMENSION Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *    COVAR(NCVM,NCVM),BETA(MMAX),P(MAXM1,MAXN)
      DATA ZERO/0D0/
      KK=MFIT+1
      DO 12 J=1,MA
       IHIT=0
       DO 11 K=1,MFIT
        IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
       IF (IHIT.EQ.0) THEN
        LISTA(KK)=J
        KK=KK+1
       ELSE IF (IHIT.GT.1) THEN
        STOP 'Improper set in LISTA'
       ENDIF
12     CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA'
      DO 14 J=1,MFIT
       DO 13 K=1,MFIT
        COVAR(J,K)=ZERO
13      CONTINUE
       BETA(J)=ZERO
14    CONTINUE
      DO 18 I=1,NDATA
       YM=Y(I)
       IF(MFIT.LT.MA) THEN
        DO 15 J=MFIT+1,MA
         YM=YM-A(LISTA(J))*P(LISTA(J),I)
15       CONTINUE
       ENDIF
       SIG2I=1./SIG(I)**2
       DO 17 J=1,MFIT
        WT=P(LISTA(J),I)*SIG2I
        DO 16 K=1,J
         COVAR(J,K)=COVAR(J,K)+WT*P(LISTA(K),I)
16       CONTINUE
        BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18     CONTINUE
      IF (MFIT.GT.1) THEN
       DO 21 J=2,MFIT
        DO 19 K=1,J-1
         COVAR(K,J)=COVAR(J,K)
19       CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
       A(LISTA(J))=BETA(J)
22     CONTINUE
      CHISQ=ZERO
      DO 24 I=1,NDATA
       SUM=ZERO
       DO 23 J=1,MA
        SUM=SUM+A(J)*P(J,I)
23      CONTINUE
       CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24     CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=1000,MAXM=10,MAXM1=MAXM+1)
      PARAMETER (NMAX=MAXM1)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DATA ZERO,ONE/0D0,1D0/
      DO 11 J=1,N
       IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
       BIG=ZERO
       DO 13 J=1,N
        IF(IPIV(J).NE.1)THEN
         DO 12 K=1,N
          IF (IPIV(K).EQ.0) THEN
           IF (ABS(A(J,K)).GE.BIG)THEN
            BIG=ABS(A(J,K))
            IROW=J
            ICOL=K
           ENDIF
          ELSE IF (IPIV(K).GT.1) THEN
           STOP 'Singular matrix'
          ENDIF
12        CONTINUE
        ENDIF
13      CONTINUE
       IPIV(ICOL)=IPIV(ICOL)+1
       IF (IROW.NE.ICOL) THEN
        DO 14 L=1,N
         DUM=A(IROW,L)
         A(IROW,L)=A(ICOL,L)
         A(ICOL,L)=DUM
14       CONTINUE
        DO 15 L=1,M
         DUM=B(IROW,L)
         B(IROW,L)=B(ICOL,L)
         B(ICOL,L)=DUM
15       CONTINUE
       ENDIF
       INDXR(I)=IROW
       INDXC(I)=ICOL
       IF (A(ICOL,ICOL).EQ.ZERO) STOP 'Singular matrix.'
       PIVINV=ONE/A(ICOL,ICOL)
       A(ICOL,ICOL)=ONE
       DO 16 L=1,N
        A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
       DO 17 L=1,M
        B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
       DO 21 LL=1,N
        IF(LL.NE.ICOL)THEN
         DUM=A(LL,ICOL)
         A(LL,ICOL)=ZERO
         DO 18 L=1,N
          A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18        CONTINUE
         DO 19 L=1,M
          B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19        CONTINUE
        ENDIF
21      CONTINUE
22     CONTINUE
      DO 24 L=N,1,-1
       IF(INDXR(L).NE.INDXC(L))THEN
        DO 23 K=1,N
         DUM=A(K,INDXR(L))
         A(K,INDXR(L))=A(K,INDXC(L))
         A(K,INDXC(L))=DUM
23       CONTINUE
       ENDIF
24    CONTINUE
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DATA ZERO/0D0/
      DO 12 J=1,MA-1
       DO 11 I=J+1,MA
        COVAR(I,J)=ZERO
11      CONTINUE
12     CONTINUE
      DO 14 I=1,MFIT-1
       DO 13 J=I+1,MFIT
        IF(LISTA(J).GT.LISTA(I)) THEN
         COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
        ELSE
         COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
        ENDIF
13      CONTINUE
14     CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
       COVAR(1,J)=COVAR(J,J)
       COVAR(J,J)=ZERO
15     CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
       COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16     CONTINUE
      DO 18 J=2,MA
       DO 17 I=1,J-1
        COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18     CONTINUE
      RETURN
      END
