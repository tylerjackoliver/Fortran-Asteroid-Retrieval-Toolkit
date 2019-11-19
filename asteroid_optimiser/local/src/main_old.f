C
C     Test problem for BOBYQA, the objective function being the sum of
C     the reciprocals of all pairwise distances between the points P_I,
C     I=1,2,...,M in two dimensions, where M=N/2 and where the components
C     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
C     defines the M points P_I. The initial X gives equally spaced points
C     on a circle. Four different choices of the pairs (N,NPT) are tried,
C     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
C     minimum that is not global occurs in both the N=10 cases. The details
C     of the results are highly sensitive to computer rounding errors. The
C     choice IPRINT=2 provides the current X and optimal F so far whenever
C     RHO is reduced. The bound constraints of the problem require every
C     component of X to be in the interval [-1,1].
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER N = 2
      DIMENSION X(N),XL(N),XU(N),W(500000)
      IPRINT=2
      MAXFUN=500000
      RHOBEG=86400.D0 * 3.D0
      RHOEND=3600.D0
      M=5
   10 N=2*M
      XL(1) = 
      CALL BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   50 CONTINUE
      M=M+M
      IF (M .LE. 10) GOTO 10
      STOP
      END
