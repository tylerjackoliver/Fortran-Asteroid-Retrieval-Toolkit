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

      DIMENSION X(2),XL(2),XU(2),W(500000)

      INTEGER I,J,K

C     OPTIMISER-SPECIFIC VALUES 

      TRANSFER_EPOCH = 0.0
      N=2
      IPRINT=2
      MAXFUN=500000

C     INITIALISE MAXIMUM AND MINIMUM STEP-SIZES FOR FINDING LOCAL MINIMA

      RHOBEG=86400.D0 * 100.D0
      RHOEND=3600.D-5

C     N + 2 < NPT < (N+2)(N+1)/2
      NPT=5

C     INITIALISE TRANSFER EPOCH 
      call FURNSH('naif0008.tls')
      call STR2ET('22 Sep 2028 00:00', TRANSFER_EPOCH)
      call UNLOAD('naif0008.tls')

C     INITIALISE BOUNDS
C     XL = [TRANSFER_EPOCH, TRANSFER_TIME]_LOWER
c     XU = [TRANSFER_EPOCH, TRANSFER_TIME]_UPPER

      XL(1)=TRANSFER_EPOCH*.99D0
      XL(2)=1241.d0*86400.d0*0.9D0
      XU(1)=TRANSFER_EPOCH*1.01
      XU(2)=1241.d0*86400.d0*1.1

C     SET INITIAL GUESS

      ! X(1)=TRANSFER_EPOCH*1.00001d0                       ! STILL CONVERGES WHEN THROWN OFF
      ! X(2)=1241.D0*86400*1.1D0                            ! STILL CONVERGES WHEN THROWN OFF
      X(1)=TRANSFER_EPOCH
      X(2)=1241.D0 * 86400.D0
      CALL BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W)

      END
