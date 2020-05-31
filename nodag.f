      SUBROUTINE NODAG(N, SIGMA, A, LAMBD, EPS, ALPHA, MAXITR)
c     Find a sparse parametrization of the inverse covariance
c     matrix as  A A**t using a proximal gradient algorithm 
c     for the problem
c               minimize   -LOGLIK(AA**t | SIGMA) + LAMBD * ||A||_1 
c     gherardo varando (2020) <gherardo.varando[at]gmail.com>
      INTEGER N,MAXITR
      DOUBLE PRECISION SIGMA(N,N),A(N,N),LAMBD,EPS,ALPHA      
ccccc   f2py  
cf2py intent(in)  n
cf2py intent(in)  sigma
cf2py intent(out) a 
cf2py double precision intent(in)  lambd
cf2py double precision intent(in,out)  :: eps = 0.0001
cf2py double precision intent(in,out)  :: alpha = 0.5
cf2py integer intent(in,out)  :: maxitr = 100 
ccccc
c     on entry
c        N      integer
c               size of the problem 
c        SIGMA  double precision (N,N)
c               empirical covariance matrix 
c        A      double precision (N,N) 
c        LAMBD double precision 
c               penalization coefficient 
c        EPS    double precision 
c               tollerance for termination
c        ALPHA  double precision 
c               coefficient for line search 
c        MAXITR integer
c               maximum number of iterations 
c     on exit
c        A      results of the optimization 
c        EPS    last difference in objective function   
c        ALPHA  minus log-likelhood on exit
c        MAXITR number of iterations 
c 
c     internal variables
      INTEGER I,J,ITR, IPIV(N),INFO
      DOUBLE PRECISION F,FNW,TMP(N,N),GRD(N,N), D(N,N), WORK(2*N),
     *                 ONE, ZERO, STEP, G, GNW, AOLD(N,N), DIFF
      ITR = 0
      ONE = 1.0
      ZERO = 0.0
c     copy A in TMP
      DO 20 J = 1,N 
         DO 10 I = 1,N
            A(I,J) = 0
            TMP(I,J) = 0
  10     CONTINUE          
         A(J,J) = 1
         TMP(J,J) = 1
  20  CONTINUE 
c     compute initial objective function
c     compute TMP = PLU 
      CALL DGETRF(N, N, TMP, N, IPIV,INFO)  
      IF (INFO .GT.0) GOTO 900
c     compute D = SIGMA * A
      CALL DSYMM("L","U", N, N, ONE, SIGMA, N, A, N, ZERO, D, N)  
      F = 0
      G = 0
c     compute F = -2*LOG(DET(A)) + trace(A**t SIGMA A) 
c             and G = ||A||_1
      DO 30 J=1,N
         F = F -  2 * LOG(ABS(TMP(J,J)))  
         DO 25 I=1,N
         F = F + A(I,J)*D(I,J)
         IF (I .NE. J) G = G + LAMBD * ABS(A(I,J))
  25     CONTINUE
  30  CONTINUE
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute TMP = TMP**-1
      CALL DGETRI(N, TMP, N, IPIV, WORK, 2*N, INFO)
c     compute GRD = 2*(D - TMP**t)
      DO 35 J=1,N
         DO 32 I=1,N
            GRD(I,J) = 2*(D(I,J) - TMP(J,I))
  32     CONTINUE   
  35  CONTINUE
c     copy old A before starting line search 
      DO 90 J = 1,N 
         DO 80 I = 1,N
             AOLD(I,J) = A(I,J)
  80     CONTINUE          
  90  CONTINUE 
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N 
         DO 100 I = 1,N
            A(I,J) = AOLD(I,J) - STEP * GRD(I,J) 
  100    CONTINUE
            IF (A(J,J) .LE. 0) THEN
                   STEP = STEP * ALPHA
                   GOTO 600 
            ENDIF
  110 CONTINUE
c     soft thresholding
      DO 130 J =1,N
         DO 120 I=1,N
            IF (I .NE. J) THEN
                A(I,J) = SIGN(ONE,A(I,J))*(ABS(A(I,J))-STEP*LAMBD) 
                IF (ABS(A(I,J)) .LE. STEP*LAMBD) THEN
                         A(I,J) = 0
                ENDIF
            ENDIF
            TMP(I,J) = A(I,J)
 120     CONTINUE
 130  CONTINUE
c     compute FNW, objective function in new A
c     compute TMP = PLU 
      CALL DGETRF(N, N, TMP, N, IPIV,INFO)  
      IF (INFO .GT. 0) THEN
          STEP = ALPHA * STEP
          GOTO 600
      ENDIF
c     compute D = SIGMA * A
      CALL DSYMM("L","U", N, N, ONE, SIGMA, N, A, N, ZERO, D, N)  
      FNW = 0
      GNW = 0
      DIFF = 0
c     compute F = -2*LOG(DET(A)) + trace(A**t SIGMA A) 
c             and G = ||A||_1
      DO 150 J=1,N
         FNW = FNW -  2 * LOG(ABS(TMP(J,J)))  
         DO 140 I=1,N
            DIFF = DIFF + ((A(I,J) - AOLD(I,J))**2)/(2*STEP) + 
     *             (A(I,J) - AOLD(I,J)) * GRD(I,J)  
            FNW = FNW + A(I,J)*D(I,J)
            IF (I .NE. J) GNW = GNW + LAMBD * ABS(A(I,J))
  140    CONTINUE
  150  CONTINUE
c     line search condition  
      IF (FNW .GT. F + DIFF .OR. FNW  + GNW .GT. F + G ) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F+G-FNW-GNW).LE.EPS).OR.(ITR.GE.MAXITR))THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F + G - FNW - GNW)  
         MAXITR = ITR
         GOTO 900 
      ENDIF  
c     update value of objective function and repeat
      F = FNW
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of NODAG
      END
c
c
