      SUBROUTINE NODAG(N, SIGMA, A, LAMBD, EPS, ALPHA, MAXITR)
c     NODAG routine (version 0.0.2)  
c     gherardo varando (2020) <gherardo.varando[at]gmail.com>
c
c     Find a sparse parametrization of the inverse covariance
c     matrix as  A A**t using a proximal gradient algorithm 
c     for the problem
c             minimize   -LOGLIK(AA**t | SIGMA) + LAMBD * ||A||_1 
c             A invertible
      IMPLICIT NONE
c     integer variables 
      INTEGER N,MAXITR
c     double precision variables
      DOUBLE PRECISION SIGMA(N,N),A(N,N),LAMBD,EPS,ALPHA      
ccccc  comment block f2py  
cf2py intent(in)  n
cf2py intent(in)  sigma
cf2py intent(out) a 
cf2py double precision intent(in)  lambd
cf2py double precision intent(in,out)  :: eps = 0.0001
cf2py double precision intent(in,out)  :: alpha = 0.5
cf2py integer intent(in,out)  :: maxitr = 100 
ccccc end comment block f2py
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
c     
c     external subroutine from LAPACK 
      EXTERNAL DGETRI, DGETRF, DSYMM
c     intrinsic functions
      INTRINSIC ABS, LOG, SIGN
c     internal variables
      INTEGER I,J,ITR, IPIV(N),INFO
      DOUBLE PRECISION F,FNW,TMP(N,N),GRD(N,N), D(N,N), WORK(2*N),
     *                 ONE, ZERO, STEP, G, GNW, AOLD(N,N), DIFF
      PARAMETER ( ONE = 1.0E+0 )
      PARAMETER ( ZERO = 0.0E+0 )
      ITR = 0
      F = 0.0E+0
c     initialize penalty for A = identity 
      G = LAMBD * N
c     initialize matrices and F = trace(SIGMA)  
      DO 20 J = 1,N 
         DO 10 I = 1,N
            A(I,J) = 0.0E+0
            TMP(I,J) = 0.0E+0
            D(I,J) = SIGMA(I,J)
  10     CONTINUE          
         A(J,J) = 1.0E+0
         TMP(J,J) = 1.0E+0
         IPIV(J) = J
         F = F + SIGMA(J,J)
  20  CONTINUE 
c     compute initial objective function in A = identity
c     compute TMP = PLU 
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute TMP = TMP**-1
      CALL DGETRI(N, TMP, N, IPIV, WORK, 2*N, INFO)
c     compute GRD = 2*(D - TMP**t)
c     and copy A before starting line search 
      DO 40 J=1,N
         DO 30 I=1,N
            GRD(I,J) = 2*(D(I,J) - TMP(J,I))
            AOLD(I,J) = A(I,J)
  30     CONTINUE   
  40  CONTINUE
      STEP = 1
c     line search loop
  600 CONTINUE     
c     gradient step
c     and soft thresholding
      DO 110 J =1,N
         DO 100 I=1,N
            A(I,J) = AOLD(I,J) - STEP * GRD(I,J) 
            A(I,J) = SIGN(ONE,A(I,J))*(ABS(A(I,J))-STEP*LAMBD) 
            IF (ABS(A(I,J)) .LE. STEP*LAMBD) THEN
               A(I,J) = 0.0E+0
            ENDIF
            TMP(I,J) = A(I,J)
 100     CONTINUE
 110  CONTINUE
c     compute TMP = PLU 
      CALL DGETRF(N, N, TMP, N, IPIV,INFO)  
      IF (INFO .GT. 0) THEN
          STEP = ALPHA * STEP
          GOTO 600
      ENDIF
c     compute D = SIGMA * A
      CALL DSYMM("L","U", N, N, ONE, SIGMA, N, A, N, ZERO, D, N)  
c     compute FNW, GNW in new A
      FNW = 0.0E+0
      GNW = 0.0E+0
      DIFF = 0.0E+0
c     compute F = -2*LOG(DET(A)) + trace(A**t SIGMA A) 
c             and G = ||A||_1
      DO 150 J=1,N
         FNW = FNW -  2 * LOG(ABS(TMP(J,J)))  
         DO 140 I=1,N
            DIFF = DIFF + ((A(I,J) - AOLD(I,J))**2)/(2*STEP) + 
     *             (A(I,J) - AOLD(I,J)) * GRD(I,J)  
            FNW = FNW + A(I,J)*D(I,J)
            GNW = GNW + LAMBD * ABS(A(I,J))
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
         EPS = F + G - FNW - GNW 
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
