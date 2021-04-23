      PROGRAM POLYFIT
      
C     Reads file 'coeff.dat' and generates E(v),P(v), and B(V) for N=3
C     E(v)=Sum[a(n) v^(-2n/3),{n,0,N}]

      DOUBLE PRECISION X,E,E1,E2,A(3),P,BM, ESTART
      DOUBLE PRECISION VOL,VOL23,TWO3,CONST,LAT
      INTEGER STRUC,ISVOL
      DATA CONST/1.47105164D2/,TWO3/6.666666666666666667D-1/
      PARAMETER (N=3)

C     STRUCT=  1-bcc, 2-fcc, 3-sc
C     ISVOL(volume units?)= 1-yes, 0-no  
      OPEN (1,FILE='coeff.inp',BLANK='ZERO',STATUS='UNKNOWN')
      READ(1,*) STRUC,ISVOL
      READ(1,*) (A(I),I=0,N)
      READ(1,*) XMIN,XMAX,DELX 
	ESTART=A(0)
      CLOSE(1)
   
      OPEN(UNIT=2,FILE='energyfit.out',STATUS='NEW')
      OPEN(UNIT=3,FILE='bulkmodfit.out',STATUS='NEW')
      OPEN(UNIT=4,FILE='pressfit.out',STATUS='NEW')
      
      X=XMIN
      DO WHILE(X.LE.XMAX)
         IF(ISVOL.EQ.1) THEN 
            VOL=X
            LAT=VOL**(1D0/3D0)
            IF(STRUC.EQ.1) THEN 
               LAT=LAT*2D0
            ELSE IF(STRUC.EQ.2) THEN 
               LAT=LAT*4D0
            ENDIF
         ELSE IF(ISVOL.EQ.0) THEN
            VOL=X**(3D0)
            LAT=X
            IF(STRUC.EQ.1) THEN 
               VOL=VOL/2D0
            ELSE IF(STRUC.EQ.2) THEN 
               VOL=VOL/4D0
            ENDIF
         ENDIF
         
         VOL23=VOL**(-TWO3)
c     E(V) in Rydberg
         E =ESTART+A(1)*VOL23+A(2)*VOL23*VOL23+A(3)*VOL23*VOL23*VOL23
C     E'(V)
         E1=A(1)+2.0*VOL23*A(2)+3.0*VOL23*VOL23*A(3)
C     E''(V)
         E2=2.0*A(2)+6.0*VOL23*A(3)
         BM=CONST*(VOL23**(2.5D0))*(4.0*VOL23*E2+10D0*E1)/9D0
c     P(V) in Mbar
         P=(TWO3*VOL23*E1/VOL)*CONST

c     Write out Energy
         if (ISVOL.eq.0) WRITE(2,38) LAT,VOL,E
         if (ISVOL.eq.1) WRITE(2,48) VOL,E
 38      FORMAT(1X,3F18.8)
 48      FORMAT(1X,2F18.8)
c     Write out Bulk Modulus
         if (ISVOL.eq.0) WRITE(3,38) LAT,VOL,BM
         if (ISVOL.eq.1) WRITE(3,48) VOL,BM
c     Write out Pressure
         if (ISVOL.eq.0) WRITE(4,38) LAT,VOL,P
         if (ISVOL.eq.1) WRITE(4,48) VOL,P
         
         X=X+DELX
      ENDDO
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      STOP
      END
