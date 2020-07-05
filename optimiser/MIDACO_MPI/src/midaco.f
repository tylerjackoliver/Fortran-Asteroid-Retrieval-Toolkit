CCCCCCCCCCCCCCCCCCCCCCCC MIDACO FORTRAN HEADER CCCCCCCCCCCCCCCCCCCCCCCCC
C                           
C     _|      _|  _|_|_|  _|_|_|      _|_|      _|_|_|    _|_|    
C     _|_|  _|_|    _|    _|    _|  _|    _|  _|        _|    _|  
C     _|  _|  _|    _|    _|    _|  _|_|_|_|  _|        _|    _|  
C     _|      _|    _|    _|    _|  _|    _|  _|        _|    _|  
C     _|      _|  _|_|_|  _|_|_|    _|    _|    _|_|_|    _|_|  
C
C                                          Version 6.0 (Limited)
C                                                           
C    MIDACO - Mixed Integer Distributed Ant Colony Optimization
C    ----------------------------------------------------------
C
C    MIDACO is a general solver for single- and multi-objective 
C    Mixed Integer Non-Linear Programs (MINLP's) of the form:
C
C
C      Minimize     F_1(X),... F_O(X)  where X(1,...N-NI)   is CONTINUOUS
C                                      and   X(N-NI+1,...N) is DISCRETE
C
C      subject to   G_j(X)  =  0   (j=1,...ME)      equality constraints
C                   G_j(X) >=  0   (j=ME+1,...M)  inequality constraints
C
C      and bounds   XL <= X <= XU
C
C
C    MIDACO is a (heuristic) global optimization solver that approximates 
C    a solution 'X' to the above displayed optimization problem. MIDACO 
C    is based on an extended Ant Colony Optimization framework (see [1]) in 
C    combination with the Oracle Penalty Method (see [2]) for constraints 'G(X)'
C    and a backtracking line search algorithm for local refinements.
C
C    MIDACO is a derivative free black box solver and can handle any kind
C    of non-smooth and highly nonlinear functions. MIDACO can also handle 
C    functions with moderate stochastic noise. For integer and mixed-integer 
C    problems, MIDACO evaluates integer variables only at true integer points. 
C    Thus, MIDACO does not perform a relaxation or other surrogate technique.
C
C    MIDACO does not require any user specified parameter tuning as it can 
C    run completely on 'Autopilot' (all parameter set equal to zero). 
C    Optionally, the user can adjust the MIDACO performance by setting 
C    some parameters explained below.
C
C    In case of mixed integer problems, the continuous variables are stored 
C    first in 'X(1,...N-NI)', while the discrete (also called integer or 
C    combinatorial) variables are stored last in 'X(N-NI+1,...N)'.
C    As an example consider:
C
C       X = (  0.1234,  5.6789,  1.0000,  2.0000,  3.0000)   
C
C       where 'N' = 5 and 'NI' = 3
C
C    Note that all 'X' is of type double precision. Equality and inequality 
C    constraints are handled in a similar way. The vector 'G' stores at first 
C    the 'ME' equality constraints and behind those, the remaining 'M-ME' 
C    inequality constraints are stored. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C    List of MIDACO subroutine arguments
C    -----------------------------------
C
C    P  :   Parallelization Factor. If no parallelization is desired, set P = 1.
C
C    O  :   Number of objective functions.
C
C    N  :   Number of optimization variables in total (continuous and integer ones). 
C           'N' is the dimension of the iterate 'X' with X = (X_1,...,X_N).
C
C    NI :   Number of integer optimization variables. 'NI' <= 'N'.
C           Integer (discrete) variables must be stored at the end of 'X'.
C     
C    M  :   Number of constraints in total (equality and inequality ones).
C           'M' is the dimension of a constraint vector 'G' with G = (G_1,...,G_M).
C
C    ME :   Number of equality constraints. 'ME' <= 'M'.
C           Equality constraints are stored in the beginning of 'G'. 
C           Inequality constraints are stored in the end of 'G'.
c
C    X(P*N) :  Array containing the iterate 'X'. For P=1 (no parallelization)
C              'X' stores only one iterate and has length 'N'. For P>1 
C              'X' contains several iterates, which are stored one after
C              another.
C
C    F(P*O) :  Array containing the objective function values 'F' corresponding
C              to the iterates 'X'. For P=1 (no parallelization), 'F' is has lenth 'O'.
C              For P>1 F has length 'P*O' and stores the vectors of objectives 
C              corresponding to 'X' one after another.
C
C    G(P*M) :  Array containing the constraint values 'G'.For P=1 (no parallelization) 
C              'G' has length 'M'. For P>1 'G' has length 'P*M' and stores the vectors
C              of constraints corresponding to 'X' one after another.
C
C    XL(N) :   Array containing the lower bounds for the iterates 'X'.
C              Note that for integer dimensions (i > N-NI) the bounds should also be 
C              discrete, e.g. X(i) = 1.0000.
C
C    XU(N) :   Array containing the upper bounds for the iterates 'X'.
C              Note that for integer dimensions (i > N-NI) the bounds should also be 
C              discrete, e.g. X(i) = 1.0000.
C
C    IFLAG :   Information flag used by MIDACO. Initially MIDACO must be called with IFLAG=0.
C              If MIDACO works correctly, IFLAG values lower than 0 are used for internal 
C              communication. If MIDACO stops (either by submitting ISTOP=1 or automatically
C              by the FSTOP or ALGOSTOP parameter), an IFLAG value between 1 and 9 is returned 
C              as final message. If MIDACO detects at start-up some critical problem setup, a 
C              ***WARNING*** message is returned by IFLAG as value between 10 and 99. If
C              MIDACO detects an ***MIDACO INPUT ERROR***, an IFLAG value between 100 and 999 
C              is returned and MIDACO stops. The individual IFLAG flags are listed below.
C
C    ISTOP :   Communication flag to stop MIDACO. If MIDACO is called with 
C              ISTOP = 1, MIDACO returns the best found solution in 'X' with 
C              corresponding 'F' and 'G'. 
C
C    PARAM() :  Array containing 13 parameters that can be selected by the user to adjust MIDACO. 
C               (See the user manual for a more detailed description of individual parameters) 
C
C    RW(LRW) :  Real workarray (Type: double precision) of length 'LRW'
C    LRW :      Length of 'RW'. 'LRW' must be greater or equal to:
C               
C                    120*N+20*M+20*O+20*P+P*(M+2*O)+O*O+5000
C
C    IW(LIW) :  Integer workarray (Type: long integer) of length 'LIW'
C    LIW :      Length of 'IW'. 'LIW' must be greater or equal to:
C
C                    3*N+P+1000
C
C    PF(LPF) :  Real array, containing the approximation of the pareto front.
C    LPF :      Length of 'PF'. 'LPF' must be greater or equal to:
C
C                     (O+M+N)*PARETOMAX + 1
C 
C               where PARETOMAX = 1000 by default. In case PARAM(10) is not equal zero,
C               PARETOMAX equals the absolute value of PARAM(10).
C
C    KEY :  License-Key for MIDACO. Note that any licensed copy of MIDACO comes with an 
C           individual 'KEY' determining the license owner and its affiliation.  
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C    List of IFLAG messages
C    ----------------------
C
C     Final Messages:
C     ---------------
C     IFLAG = 1 : Feasible solution,   MIDACO was stopped by MAXEVAL or MAXTIME
C     IFLAG = 2 : Infeasible solution, MIDACO was stopped by MAXEVAL or MAXTIME
C     IFLAG = 3 : Feasible solution,   MIDACO stopped automatically by ALGOSTOP
C     IFLAG = 4 : Infeasible solution, MIDACO stopped automatically by ALGOSTOP
C     IFLAG = 5 : Feasible solution,   MIDACO stopped automatically by EVALSTOP
C     IFLAG = 6 : Infeasible solution, MIDACO stopped automatically by EVALSTOP
C     IFLAG = 7 : Feasible solution,   MIDACO stopped automatically by FSTOP
C
C     WARNING - Flags:
C     ----------------
C     IFLAG = 51 : Some X(i)  is greater/lower than +/- 1.0D+16 (try to avoid huge values)
C     IFLAG = 52 : Some XL(i) is greater/lower than +/- 1.0D+16 (try to avoid huge values)
C     IFLAG = 53 : Some XU(i) is greater/lower than +/- 1.0D+16 (try to avoid huge values)
C     IFLAG = 71 : Some XL(i) = XU(i) (fixed variable)
C     IFLAG = 81 : F(1) has value NaN for starting point X
C     IFLAG = 82 : Some G(X) has value NaN for starting point X 
C     IFLAG = 91 : FSTOP is greater/lower than +/- 1.0D+16
C     IFLAG = 92 : ORACLE is greater/lower than +/- 1.0D+16
C
C     ERROR - Flags:
C     --------------
C     IFLAG = 100 :   P   <= 0   or   P  > 1.0D+99
C     IFLAG = 101 :   O   <= 0   or   O  > 1000000000
C     IFLAG = 102 :   N   <= 0   or   N  > 1.0D+99
C     IFLAG = 103 :   NI  <  0
C     IFLAG = 104 :   NI  >  N
C     IFLAG = 105 :   M   <  0   or   M  > 1.0D+99
C     IFLAG = 106 :   ME  <  0
C     IFLAG = 107 :   ME  >  M
C     IFLAG = 201 :   some X(i)  has type NaN
C     IFLAG = 202 :   some XL(i) has type NaN
C     IFLAG = 203 :   some XU(i) has type NaN
C     IFLAG = 204 :   some X(i) < XL(i)
C     IFLAG = 205 :   some X(i) > XU(i)
C     IFLAG = 206 :   some XL(i) > XU(i)
C     IFLAG = 301 :   PARAM(1) < 0   or   PARAM(1) > 1.0D+99
C     IFLAG = 302 :   PARAM(2) < 0   or   PARAM(2) > 1.0D+99
C     IFLAG = 303 :   PARAM(3) greater/lower than +/- 1.0D+99
C     IFLAG = 304 :   PARAM(4) < 0   or   PARAM(4) > 1.0D+99
C     IFLAG = 305 :   PARAM(5) greater/lower than +/- 1.0D+99
C     IFLAG = 306 :   PARAM(6) not a discrete value or   PARAM(6) > 1.0D+99
C     IFLAG = 307 :   PARAM(7) < 0   or   PARAM(7) > 1.0D+99
C     IFLAG = 308 :   PARAM(8) < 0   or   PARAM(8) > 100
C     IFLAG = 309 :   PARAM(7) < PARAM(8)
C     IFLAG = 310 :   PARAM(7) > 0 but PARAM(8) = 0
C     IFLAG = 311 :   PARAM(8) > 0 but PARAM(7) = 0
C     IFLAG = 312 :   PARAM(9) greater/lower than +/- 1.0D+99
C     IFLAG = 321 :   PARAM(10) >= 1.0D+99
C     IFLAG = 322 :   PARAM(10) not a discrete value
C     IFLAG = 331 :   PARAM(11) < 0 or > 0.5
C     IFLAG = 344 :   LPF too small. LPF must be at least (O+M+N)*PARETOMAX + 1,
C                     where PARETOMAX = 1000 by default. See also PARAM(10).
C     IFLAG = 347 :   PARAM(5) > 0 but PARAM(5) < 1
C     IFLAG = 348 :   PARAM(5): Optional EVALSTOP precision appendix > 0.5
C     IFLAG = 350 :   PARAM(12) > 1 or < 1 but not a discrete value
C     IFLAG = 351 :   PARAM(13) < 0   or   PARAM(13) > 3
C     IFLAG = 352 :   PARAM(13) not a discrete value
C     IFLAG = 399 :   Some PARAM(i) has type NaN
C     IFLAG = 401 :   ISTOP < 0 or ISTOP > 1
C     IFLAG = 402 :   Starting point does not satisfy all-different constraint
C     IFLAG = 501 :   Double precision work space size LRW is too small.  
C                     Increase size of RW array. RW must be at least of 
C                     size LRW = 120*N+20*M+20*O+20*P+P*(M+2*O)+O*O+5000
C
C     IFLAG = 502 :   Internal LRW check error
C     IFLAG = 601 :   Integer work space size LIW is too small.
C                     Increase size of IW array. IW must be at least of 
C                     size LIW = 3*N+P+1000
C
C     IFLAG = 602 :   Internal LIW check error 
C     IFLAG = 701 :   Input check failed. MIDACO must be called initially with IFLAG = 0
C     IFLAG = 801 :   P > PMAX (user must increase PMAX in the MIDACO source code) 
C     IFLAG = 881 :   Integer part of X contains continues (non discrete) values
C     IFLAG = 882 :   Integer part of XL contains continues (non discrete) values
C     IFLAG = 883 :   Integer part of XU contains continues (non discrete) values
C     IFLAG = 900 :   Invalid or corrupted LICENSE-KEY
C     IFLAG = 999 :   N > 4. The free test version is limited up to 4 variables.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C    References:
C    -----------
C
C    [1] Schlueter, M., Egea, J. A., Banga, J. R.: 
C        "Extended ant colony optimization for non-convex mixed integer nonlinear programming", 
C        Computers and Operations Research (Elsevier), Vol. 36 , Issue 7, Page 2217-2229, 2009.
C
C    [2] Schlueter M., Gerdts M.: "The oracle penalty method",
C        Journal of Global Optimization (Springer), Vol. 47, Issue 2, Page 293-325, 2010.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C    Author (C) :   Dr. Martin Schlueter
C                   Information Initiative Center,
C                   Division of Large Scale Computing Systems,
C                   Hokkaido University, JAPAN.
C
C    Email :        info@midaco-solver.com
C    URL :          http://www.midaco-solver.com
C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE MIDACO( P, O, N, NI, M, ME, X, F, G, XL, XU, IFLAG, 
     &                   ISTOP, PARAM, RW, LRW, IW, LIW, PF, LPF, KEY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER P, O, N, NI, M, ME, IFLAG, LRW, IW, LIW, LPF, ISTOP
      DOUBLE PRECISION X, F, G, XL, XU, PARAM(13), RW, PF
      CHARACTER KEY*60
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION X(P*N),F(P*O),G(P*M),XL(N),XU(N),RW(LRW),IW(LIW),PF(LPF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER eA,eB,eDG,eDF,eU,eN,eIE,eM,eFM, OFFSET
      DATA    eA,eB,eDG,eDF,eU,eN,eIE,eM,eFM /0,0,0,0,0,0,0,0,0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Check number of variables for the limited version. Note that 
C     removing this check will cause MIDACO to not work properly on 
C     problems with N > 4.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N.GT.4)THEN
        IFLAG = 999
        ISTOP = 1
        RETURN
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      IF( IFLAG .EQ. 0 )THEN

         CALL PRECHECK(P,O,N,M,LRW,LIW,LPF,PF,PARAM,IFLAG,ISTOP)

         OFFSET = P+O+O+N+M+200

C        Setting up RW workspace entry pointers
         eA  = OFFSET + 107*n+6*m+7*o+5*p+616        
         eB  = OFFSET + eA  + P 
         eDG = OFFSET + eB  + P
         eDF = OFFSET + eDG + P*(M+2*O)+1 
         eU  = OFFSET + eDF + P
         eN  = OFFSET + eU  + O
         eIE = OFFSET + eN  + O
         eM  = OFFSET + eIE + P
         efM = OFFSET + eM  + P         
         
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL MIDACO_KERNEL_DRIVER(P,O,N,NI,M,ME,F,G,X,XL,XU, 
     &                          IFLAG,ISTOP,PARAM,RW,LRW,IW,LIW,PF,
     &                          LPF,eA,eB,eDG,eDF,eU,eN,eIE,eM,efM,KEY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
              subroutine i409(m,i8,g,i5,i2,i4,i6,i43)         
                                                   implicit none        
            integer m,i8                                                
            double precision g(*),i5(*),i2(*)                           
        integer i,j,k,i450,i10,i43, i6(*), i301,ok,t                
                      double precision o25,i4(*),z                      
                                              do k = 1,m                
                   i6(i43+k-1) = k                                    
             enddo                                                      
           do k = 1,m                                                   
                  j             = int(dble(k)*o25(i4)) + 1              
                              i6(i43+k-1) = i6(i43+j-1)             
              i6(i43+j-1) = k                                         
                                                                 enddo  
                                                   do t=1,m             
                                     i = i6(i43+t-1)                  
          if( i .le. m-i8 )goto 79                                      
                  do k=m-i8+1,m                                         
                           if( i301(g(i),g(k)).eq.1 .and. i.ne.k )then 
                      i450 = 0                                         
                i10   = 0                                               
                        do while( i450 .eq. 0 )                        
                                                     i10 = i10 + 1      
                            z = g(i) + dble(i10)                        
                            if(z.gt.i2(i)) goto 88                      
                                                            ok=0        
                                     do j=m-i8+1,m                      
                       if(i.ne.j.and.i301(z,g(j)).eq.0) ok=ok+1        
                                                         enddo          
                                             if(ok.eq.i8-1)then         
         i450 = 1                                                      
                                g(i) = z                                
                        goto 123                                        
                                                        endif           
   88                                  continue                         
                          z = g(i) - dble(i10)                          
                              if(z.lt.i5(i)) goto 99                    
                                         ok=0                           
                    do j=m-i8+1,m                                       
        if(i.ne.j.and.i301(z,g(j)).eq.0) ok=ok+1                       
                                                              enddo     
                    if(ok.eq.i8-1)then                                  
        i450 = 1                                                       
                    g(i) = z                                            
                                      goto 123                          
                                                        endif           
   99                    continue                                       
                  if( dble(i10) .gt. i2(i)-i5(i) )then                  
                                     goto 123                           
                      endif                                             
                                                               enddo    
  123                                           continue                
                                             endif                      
        enddo                                                           
   79        continue                                                   
                                                       enddo            
                                                                   end  
            subroutine i405(l,x,n,i0,i16,i48,i41,i42)   
            implicit none                                               
                       integer n,i0,i41,i42                         
          double precision l,x(*),i48(*),i16                          
                                     integer k2,i301                  
           double precision i56,i75,i17,k1                       
                                         data k2 /0/                   
            data i56,i75,k1,i17 /0.0d0,0.0d0,0.0d0,0.0d0/        
                                                 if( i42 .eq. 0 )then 
           k2 = 0                                                      
            call i403( x,n,i0,i16,i17 )                    
                       i56   = l                                      
                                        i75 = i17                   
                k1 = 0.001d0                                           
                  if( i48(5) - dble(nint(i48(5))) .gt. 0.0d0 )then  
                          k1 = i48(5) - dble(nint(i48(5)))         
                                    endif                               
                                                     goto 9             
                                                                 endif  
                                                     if( n .le. 0 )then 
               if( l .le. i56 - dabs(i56) * k1 ) goto 7            
             else                                                       
                    call i403( x,n,i0,i16,i17 )            
          if( i75 .le. i16 )then                                    
                   if( i301( i17 , i75) .eq. 1 .and.               
     &  l .le. i56 - dabs(i56) * k1) goto 7                        
                                               else                     
           if( i17 .le. i75 - i75 * k1 ) goto 7                
                                                endif                   
                                         endif                          
                 k2 = k2 + 1                                          
            if( k2 .ge. nint(i48(5)) ) goto 8                        
                                                                goto 9  
    7                      continue                                     
                   i56   = l                                          
                                    i75 = i17                       
               k2 = 0                                                  
                                               goto 9                   
    8                        continue                                   
                                                             i41 = 1  
        if( i75 .le. i16 ) i42 = 5                                
           if( i75 .gt. i16 ) i42 = 6                             
    9          return                                                   
                          end                                           
                  subroutine i403( x,n,i0,i16,i17 )        
         implicit none                                                  
                                double precision x(*),i16,i17           
                         integer i,n,i0                                 
                                                    i17 = 0.0d0         
                                                        do i = 1,n      
                                                if(i.le.i0)then         
           if( dabs(x(i)) .gt. i16) i17 = i17 + dabs(x(i))              
                                                  else                  
                              if( x(i) .lt. -i16) i17 = i17 - x(i)      
                                    endif                               
                                                  enddo                 
       end                                                              
                       subroutine o19(f,o,m,i8,n,i0,g,l,x,i5,i2, 
     &                       i42,i41,i48,i4,i32,i6,i99,           
     &             i30,i52,i50,i100,                       
     &                                      i990)                
                                             implicit none              
                 integer f,o,m,i8,n,i0,i42,i32,i6,i99,i41           
                       double precision g,l,x,i5,i2,i48(*),i4         
                 dimension g( f*m ),l( f ),x( f*n+1 ),i5(m),i2(m)       
                          dimension i4(i32),i6(i99)                     
                         character i990*60                       
                      integer i30,i52,i50,i100             
                          integer i304, i301, i,k                    
          integer i302, i303, i436              
                  integer egtraollset                                   
                                            if( i42 .ge. 100 ) return 
                               i30 = 100                               
                     if(f.le.0.or.f.gt.1.0d99)then                      
                                                   i42 = 100          
                                             goto 701                   
                                            endif                       
                                     if(m.le.0.or.m.gt.1.0d99)then      
                                                      i42 = 102       
          goto 701                                                      
                                     endif                              
                            if(i8.lt.0)then                             
                                              i42 = 103               
                                             goto 701                   
                                        endif                           
                                        if(i8.gt.m)then                 
                                        i42 = 104                     
                                       goto 701                         
                                         endif                          
                        if(n.lt.0.or.n.gt.1.0d99)then                   
                          i42 = 105                                   
                              goto 701                                  
                      endif                                             
        if(i0.lt.0)then                                                 
                                              i42 = 106               
        goto 701                                                        
                                      endif                             
                     if(i0.gt.n)then                                    
                                      i42 = 107                       
                             goto 701                                   
                                                                  endif 
                                                              do i=1,m  
         if( i304(g(i)).eq.1 )then                                    
                                              i42 = 201               
                                                  goto 701              
                       endif                                            
                         if( i304(i5(i)).eq.1 )then                   
                                   i42 = 202                          
                       goto 701                                         
       endif                                                            
                                     if( i304(i2(i)).eq.1 )then       
             i42 = 203                                                
                                              goto 701                  
                                             endif                      
             if(g(i).lt.i5(i)-1.0d-4)then                               
                    i42 = 204                                         
           goto 701                                                     
                                      endif                             
                                      if(g(i).gt.i2(i)+1.0d-4)then      
                i42 = 205                                             
                      goto 701                                          
                        endif                                           
                                 if(i5(i).gt.i2(i)+1.0d-4)then          
                       i42 = 206                                      
          goto 701                                                      
             endif                                                      
                                   enddo                                
                  if(i48(1).lt.0.0d0.or.i48(1).gt.1.0d99)then       
                    i42 = 301                                         
            goto 701                                                    
       endif                                                            
             if(i48(2).lt.0.0d0.or.i48(2).gt.1.0d99)then            
                                  i42 = 302                           
                                                           goto 701     
                                      endif                             
                  if(i48(3).gt.1.0d99.or.i48(3).lt.-1.0d99)then     
                                                     i42 = 303        
                                                   goto 701             
        endif                                                           
             if(i48(4).lt.0.0d0.or.i48(4).gt.1.0d99)then            
                  i42 = 304                                           
                                           goto 701                     
                     endif                                              
                      if(i48(5).lt.0.0d0.or.i48(5).gt.1.0d99)then   
                                             i42 = 305                
                                    goto 701                            
                                                      endif             
                                    if(i48(6).gt.1.0d99.or.           
     &          dabs(i48(6)-dble(nint(i48(6)))) .gt. 1.0d-6)then    
                                                            i42 = 306 
                 goto 701                                               
                                             endif                      
                     if(i48(7).lt.0.0d0.or.i48(7).gt.1.0d99)then    
                        i42 = 307                                     
                                  goto 701                              
                                          endif                         
              if(i48(8).lt.0.0d0.or.i48(8).gt.i30)then             
                 i42 = 308                                            
                 goto 701                                               
                                                      endif             
         if(i48(7).gt.0.0d0.and.i48(7).lt.i48(8))then             
                                                      i42 = 309       
                                                      goto 701          
                                           endif                        
      if(i48(7).gt.0.0d0.and.(i301(i48(8),0.0d0).eq.1))then        
                 i42 = 310                                            
                                                      goto 701          
          endif                                                         
            if( (i301(i48(7),0.0d0).eq.1).and.i48(8).gt.0.0d0)then 
                                            i42 = 311                 
                      goto 701                                          
                                    endif                               
              if(i48(9).gt.1.0d99.or.i48(9).lt.-1.0d99)then         
                   i42 = 312                                          
                  goto 701                                              
                                                             endif      
                if( dabs(i48(12)) .gt. 0.0d0 )then                    
        if(  dabs(i48(12)-dble(nint(i48(12)))) .gt. 1.0d-6          
     &            .and. dabs(i48(12)) .gt. 1.0d0 )then                
                                                         i42 = 350    
                                    goto 701                            
                                          endif                         
                                                    endif               
                if(i48(13).lt.0.0d0.or.i48(13).gt.3.0d0)then        
                                         i42 = 351                    
                                             goto 701                   
                                          endif                         
           if( dabs(i48(13)-dble(nint(i48(13)))) .gt. 1.0d-6)then   
                                       i42 = 352                      
                                    goto 701                            
                                  endif                                 
                       do i=1,12                                        
                if( i304(i48(i)).eq.1 )then                         
        i42 = 399                                                     
                                           goto 701                     
                                                                 endif  
                            enddo                                       
              if( i48(5) .gt. 0.0d0 .and. i48(5) .lt. 1.0d0)then    
                                                             i42 = 347
                                                 goto 701               
              endif                                                     
                 if( i48(5) - dnint(i48(5)) .lt. 0.0d0 )then        
                                 i42 = 348                            
                                                             goto 701   
        endif                                                           
             if(i41.lt.0.or.i41.gt.1)then                           
                i42 = 401                                             
                               goto 701                                 
                                                        endif           
                  if( i301(i48(13),3.0d0) .eq. 1 )then               
                            do i=m-i8+1,m                               
                                         do k=m-i8+1,m                  
               if( i .ne. k .and. dabs(g(i)-g(k)) .le. 0.1 )then        
                      i42 = 402                                       
                 goto 701                                               
                                                                  endif 
                                   enddo                                
                            enddo                                       
                                                        endif           
                                        if(o.le.1)then                  
                                    i436  = n                         
                                                   else                 
                                      i436 = n - 2*o                  
                                            endif                       
                             egtraollset = 5*(f+o+m+i436)+100         
                i52 = 2*m+(i436+2*o)+(m+5)*i30+16 + egtraollset    
                                           i50 = 31+f+m+m             
           i302 = n*f+2*n+104*m+o*o+2*o*f+6*o+3*f+522   +2*f 
      i303 = 3*m+f+92                                        
                  if(i32.lt.i302)then                        
                i42 = 502                                             
                                   goto 701                             
       endif                                                            
                     if(i99.lt.i303)then                     
                                                          i42 = 602   
                            goto 701                                    
                                  endif                                 
                                       do i=10,i52 + 2*m+n+3          
                                     i4(i) = 0.0d0                      
                                       enddo                            
                                 do i=1,i303                 
                                i6(i) = 0                               
                          enddo                                         
                          do i = m-i8+1,m                               
           if(dabs(g(i)-dnint(g(i))).gt.1.0d-3)then                     
                                                       i42 = 881      
                                  goto 701                              
                       endif                                            
                      if(dabs(i5(i)-dnint(i5(i))).gt.1.0d-3)then        
                            i42 = 882                                 
         goto 701                                                       
            endif                                                       
                   if(dabs(i2(i)-dnint(i2(i))).gt.1.0d-3)then           
                                           i42 = 883                  
      goto 701                                                          
                                                               endif    
                                           enddo                        
            call o18(i6(i50+1),i990)                   
                               i6(i50+61) = 0                         
                           do i=1,60                                    
                            i6(i50+61) =i6(i50+61) + i6(i50+i)    
                                                 enddo                  
                              if(i6(i50+61).ne.2736)then              
                                                       i42 = 900      
                                                             goto 701   
                  endif                                                 
                                                   do i=4,7             
                                            i4(i)=dble(i6(i50+i-3))   
                                           enddo                        
       i4(8)=dble(i6(i50+61))                                         
                   i100 = 0                                     
                                                         do i=1,m       
          if(g(i).gt.1.0d16.or.g(i).lt.-1.0d16)then                     
        i42 = 51                                                      
           goto 702                                                     
                        endif                                           
                if(i5(i).gt.1.0d16.or.i5(i).lt.-1.0d16)then             
                                                    i42 = 52          
                                       goto 702                         
                                       endif                            
                     if(i2(i).gt.1.0d16.or.i2(i).lt.-1.0d16)then        
                                  i42 = 53                            
                              goto 702                                  
                                         endif                          
                      if( i301(i5(i),i2(i)).eq.1 )then                 
                 i42 = 71                                             
                    goto 702                                            
                                                              endif     
                           enddo                                        
               if( i304(l(1)).eq.1 )then                              
                                  i42 = 81                            
           goto 702                                                     
                                                       endif            
                            do i = 1,n                                  
      if( i304(x(i)).eq.1 )then                                       
                 i42 = 82                                             
                                        goto 702                        
                                                       endif            
                         enddo                                          
                    if(dabs(i48(3)).gt.1.0d16)then                    
                                    i42 = 91                          
                     goto 702                                           
                                                               endif    
                               if(dabs(i48(9)).gt.1.0d16)then         
                                          i42 = 92                    
        goto 702                                                        
              endif                                                     
                                                       return           
  701                                                          continue 
                                            i41 = 1                   
                 return                                                 
  702                                  continue                         
                           i100 = 1                             
                                                           return       
                                                            end         
      subroutine i410(o,m,n,l,x,g,pl,px,pg,i449,i452,k6,i448,
     &             i306, i307, g004,                        
     &                       g005, fx13, io10,      
     & io16, io2 )                                 
                                                       implicit none    
                                integer o,m,n,i449, i452, i448        
            double precision l(*),x(*),g(*),pl(*),px(*),pg(*),k6       
                   double precision g004,fx14,k6_i438  
       integer g005,fx13(*),io10                    
       integer i,j,k,i437,i438,i453,ii,k5,i440,i441             
                        integer fx03, fx05,fx06
      double precision z,y,k4,i439,i411,fx11,fx12  
         double precision i306(*),i307(*),fx21(1000),i305 
                           double precision rio5,fx22,k2  
                   double precision fx01,io16        
       integer fx08(1000),fx04,io2         
                   double precision fx02                      
           data i453,ii,k5 /0,0,0/                                      
                                data i411 /0.0d0/                
           io2 = 0                                        
                                  if(i449.eq.0)then                    
                        i453 = 0                                        
                                                  ii = 0                
                                           k5 = 0                       
                       i411 = 0.0d0                              
                               goto 888                                 
                               endif                                    
                if(i449.ge.2)then                                      
         do i = 1,i449                                                 
                                    fx02 = 0.0d0              
                                                           do j = 1,o   
                          z = pl(o*(i-1)+j)                             
                     fx02 = fx02 + dabs(l(j)-z)     
              enddo                                                     
      if( fx02 .le. 1.0d-6 ) return                           
                                                             enddo      
         endif                                                          
                      fx14 = k6 - k6 / dble(g005)       
                                   k6_i438 = 1.0d-12                 
           i438 = 0                                                   
             fx04 = 0                                           
                                                       do j = 1,o       
        fx21(j) = max( 1.0d-8, i307(j)-i306(j) )                     
                           fx04 = fx04 + fx13(j)    
                                                        enddo           
                                do i = 1,i449                          
                                                        i437 = 0      
                                                    fx05 = 0   
         fx03  = 1                                              
                                         do j = 1,o                     
                           if( fx13(j) .eq. 1 ) goto 123            
       z = pl(o*(i-1)+j)                                                
                                                   if( l(j) .le. z )then
              if( l(j) .le. z - fx21(j)*fx14 ) i437=i437+1     
            if( l(j) .lt. z - fx21(j)*0.01d0 ) fx05 = 1        
                                                     else               
                 if( l(j) .gt. z + fx21(j)*0.001d0 ) fx03 = 0   
                                               endif                    
  123                                                         continue  
                enddo                                                   
                  if( i437       .lt.  o-fx04 .and.           
     &                           fx03 .eq.  1 .and.             
     &    fx05 .eq. 1 )then                                    
        pl(o*(i-1)+1) = -30111979.0d0                                   
                                             i438 = i438 + 1        
                                                                   endif
                       if(i437.eq.0.and.i438.eq.0) return           
            if(i437.eq.o-fx04)then                            
                        pl(o*(i-1)+1) = -30111979.0d0                   
                   i438 = i438 + 1                                  
                                                               endif    
                                                 enddo                  
                              if( i449 .le. 10 ) goto 777              
                                           if( i438 .ge. 1 ) goto 777 
       fx06  = 0                                                
                         do i=1,o                                       
             if( l(i) .gt. i307(i)   ) fx06  = 1               
                                              enddo                     
                 if( fx06 .ge. 1 )then                          
      fx22 = 0.0d0                                                      
                                                           k2  = 0.0d0 
                       do i=1,o                                         
      if(l(i).gt.i307(i) )k2 =k2 + dabs(l(i)-i307(i))/fx21(i)       
        if(l(i).lt.i306(i))fx22=fx22+ dabs(l(i)-i306(i))/fx21(i)    
                                          enddo                         
                                   if( k2 .gt. 100.0d0 * fx22 )then    
                             if( o .le. 2 )then                         
                                    io2 = 0               
       return                                                           
                                                            else        
                             if( k2 .gt. 0.01d0 )then                  
              io2 = 0                                     
                         return                                         
         endif                                                          
                                       endif                            
                                                             endif      
                   endif                                                
  777                       continue                                    
                          if(i438.eq.0) goto 888                      
         do i = 1,i449-i438                                          
                if( dabs(pl(o*(i-1)+1)+30111979.0d0).le.k6_i438)then 
                         do j = i449-i438+1,i449                    
        if(dabs(pl(o*(j-1)+1)+30111979.0d0).gt.k6_i438)then          
                                          do k=1,o                      
                  pl(o*(i-1)+k) = pl(o*(j-1)+k)                         
                                                           enddo        
                                                     do k=1,n           
                                px(n*(i-1)+k) = px(n*(j-1)+k)           
            enddo                                                       
                         do k=1,m                                       
           pg(m*(i-1)+k) = pg(m*(j-1)+k)                                
                                       enddo                            
                                           pl(o*(j-1)+1) = -30111979.0d0
                             goto 1                                     
                                             endif                      
       enddo                                                            
                                                    endif               
    1                         continue                                  
                                                             enddo      
        i449 = i449 - i438                                          
  888      continue                                                     
       io2 = 1                                            
          if(i449.eq.i452 .and. i448.eq.1) goto 999                   
                   if(i449.eq.i452)then                                
                           k4 = i305()                      
                                                        do i=1,i449    
                                         i439 = 0.0d0                 
                                                     do j=1,o           
            z = pl(o*(i-1)+j)                                           
                  i439 = i439 + abs(l(j)-z)                         
                      enddo                                             
                                if(i439.lt.k4) k4 = i439          
                          enddo                                         
       if(i453.eq.1.and.i411.gt.k4)then                         
             goto 999                                                   
                                      endif                             
                                               i440=1                 
          i441=1                                                      
      if(i453.eq.2.and.i411.gt.k4)then                          
                                                 i440 = ii            
                                                     i441 = k5        
                            endif                                       
                               i411 = i305()          
                                                    do i=i440,i449   
                                           do k=i441,i449            
             if(i.ne.k)then                                             
                                   i439 = 0.0d0                       
                                do j=1,o                                
                                                z = pl(o*(i-1)+j)       
                                                 y = pl(o*(k-1)+j)      
                                          i439 = i439 + abs(y-z)    
                            enddo                                       
                                             if(i439.lt.k4)then      
              do j=1,o                                                  
          pl(o*(i-1)+j) = l(j)                                          
            enddo                                                       
                         do j=1,n                                       
                              px(n*(i-1)+j) = x(j)                      
                                                                 enddo  
                                      do j=1,m                          
                                          pg(m*(i-1)+j) = g(j)          
       enddo                                                            
          i453 = 2                                                      
                                                        ii = i-1        
                                                             k5 = k-1   
                                   if(ii.le.0) ii = 1                   
             if(k5.le.0) k5 = 1                                         
                               goto 999                                 
                 endif                                                  
                   if(i439 .lt. i411 ) i411 = i439    
                            endif                                       
      enddo                                                             
             enddo                                                      
                                 i453 = 1                               
                                        goto 999                        
                endif                                                   
                                                    i449 = i449 + 1   
                                             do i=1,o                   
                           pl(o*(i449-1)+i) = l(i)                     
                                                        enddo           
                                                               do i=1,n 
       px(n*(i449-1)+i) = x(i)                                         
                                       enddo                            
                                     do i=1,m                           
                                        pg(m*(i449-1)+i) = g(i)        
                                            enddo                       
  999                                         continue                  
                                                       do i=1,o         
                                             fx12 = i306(i)     
            fx11  = i307(i)                                       
                           i306(i)  =  i305()              
                           i307(i)   = -i305()              
         do k=1,i449                                                   
                                                      z = pl(o*(k-1)+i) 
                          if( z .lt. i306(i) ) i306(i) = z          
                       if( z .gt. i307(i)  )then                       
                                                 i307(i)  = z          
                                                 fx08(i) = k     
                                                           endif        
                        enddo                                           
        if( rio5( fx12, i306(i)) .gt. 1.0d-16 )then
             call o9( g004, i306(i), fx12,
     &                            o, g005, io10,        
     &             io16 )                                        
                                                endif                   
         if( rio5( fx11, i307(i)) .gt. 1.0d-16 )then 
         call o9( g004, i307(i), fx11,      
     &               o, g005, io10,                     
     &                       io16 )                              
                  endif                                                 
         enddo                                                          
                                  if( i449 .lt. 4 ) return             
                                                            do j=1,o    
                           do k=1,i449                                 
                                               k4 = 0.0d0              
          do i=1,o                                                      
              if(i.ne.j) k4 = k4 + dabs( pl(o*(fx08(j)-1)+i) - 
     &                                    pl(o*(k-1)+i) ) /             
     &                          fx21(i)                                 
                                                        enddo           
                                               if( k4 .lt. 0.03d0 )then
                   fx01 = dabs( pl(o*(fx08(j)-1)+j) -
     &                          pl(o*(k-1)+j) ) /                       
     &              fx21(j)                                             
                   if( fx01 .gt. 0.03d0 )then               
        if( fx01 / max(1.0d-12,k4) .gt. 30.0d0 )then       
                                                           k4 = 0.0d0  
                                              do i=1,o                  
                 k4 = k4 + dabs(pl(o*(fx08(j)-1)+i)-l(i))      
      enddo                                                             
       if( k4 .le. 1.0d-12 ) io2 = 0                     
            call ol004(  o,n,m,pl,px,pg,i449, fx08(j))   
                                        fx11 = i307(j)            
              i307(j) = pl(o*(k-1)+j)                                  
                      fx08(j) = k                                
        call o9( g004, i307(j), fx11,       
     &                            o, g005, io10,        
     &                  io16 )                                   
           endif                                                        
                  endif                                                 
                                             endif                      
                                                              enddo     
                                   enddo                                
                                                     end                
        subroutine  ol004( o,n,m, pl,px,pg, i449, iii )         
           implicit none                                                
                                                            integer iii 
               integer o,n,m, i449, i                                  
                                      double precision pl(*),px(*),pg(*)
              if( iii .eq. i449 )then                                  
           i449 = i449 - 1                                            
                       return                                           
                                         endif                          
                                                              do i=1,o  
                       pl(o*(iii-1)+i) = pl(o*(i449-1)+i)              
                                                        enddo           
         do i=1,n                                                       
          px(n*(iii-1)+i) = px(n*(i449-1)+i)                           
                                                                   enddo
                 do i=1,m                                               
       pg(m*(iii-1)+i) = pg(m*(i449-1)+i)                              
                    enddo                                               
                            i449 = i449 - 1                           
                                                   end                  
              subroutine o4( i67, g )            
                     implicit none                                      
                          integer i67                                
         double precision g(33)                                         
                                      if( i67 .eq. 1 )then           
            g(  1) =     26.297690237039237d0                           
                              g(  2) =      0.187494163680378d0         
                            g(  3) =      0.194021736913887d0           
                g(  4) =      0.024775540335352d0                       
                   g(  5) =   1746.883434934608658d0                    
                  g(  6) =    303.151649371507347d0                     
       g(  7) =    544.719488383699286d0                                
                                g(  8) =      0.000000000000000d0       
            g(  9) =      0.988103275217050d0                           
                               g( 10) =      0.500000000000000d0        
                              g( 11) =      0.000034825114132d0         
                   g( 12) =      0.000000000000000d0                    
                         g( 13) =      0.000000000000000d0              
                             g( 14) =      0.495363052732969d0          
                   g( 15) =     25.371780665558138d0                    
                   g( 16) =      5.150978141323659d0                    
                      g( 17) =      0.000000000000000d0                 
               g( 18) =     27.000000000000000d0                        
            g( 19) =    487.000000000000000d0                           
                                      g( 20) =     11.000000000000000d0 
                    g( 21) =    121.000000000000000d0                   
                g( 22) =      9.000000000000000d0                       
                             g( 23) =     12.000000000000000d0          
                    g( 24) =      1.000000000000000d0                   
           g( 25) =      3.000000000000000d0                            
                g( 26) =     84.000000000000000d0                       
                                   g( 27) =     47.000000000000000d0    
                 g( 28) =     14.000000000000000d0                      
             g( 29) =    505.000000000000000d0                          
                    g( 30) =     11.000000000000000d0                   
                                       g( 31) =     37.000000000000000d0
                           g( 32) =  45258.000000000000000d0            
                                g( 33) =     12.000000000000000d0       
              endif                                                     
                  if( i67 .eq. 2 )then                               
             g(  1) =       1.146763245087168d0                         
                          g(  2) =       0.026056729677485d0            
          g(  3) =       0.811482535599442d0                            
                            g(  4) =       0.940201142299118d0          
                    g(  5) =    4032.292473666903788d0                  
                                     g(  6) =     202.642029563531736d0 
          g(  7) =    1313.372702421069562d0                            
              g(  8) =       0.601346603393201d0                        
                 g(  9) =       0.383298858515245d0                     
                              g( 10) =       0.099859862564629d0        
                 g( 11) =       0.000257589716335d0                     
            g( 12) =       0.000000000000000d0                          
               g( 13) =       0.013130363400112d0                       
                               g( 14) =       0.285582894290092d0       
              g( 15) =     358.7106127644310d0                          
       g( 16) =      11.186279370841920d0                               
        g( 17) =       0.000000000000000d0                              
                g( 18) =      70.000000000000000d0                      
                  g( 19) =       3.000000000000000d0                    
                              g( 20) =      38.000000000000000d0        
                     g( 21) =      39.000000000000000d0                 
                g( 22) =      12.000000000000000d0                      
                       g( 23) =       7.000000000000000d0               
               g( 24) =       3.000000000000000d0                       
                                  g( 25) =       2.000000000000000d0    
                            g( 26) =     123.000000000000000d0          
                        g( 27) =    1271.000000000000000d0              
                                g( 28) =      25.000000000000000d0      
                               g( 29) =     594.000000000000000d0       
                                  g( 30) =       4.000000000000000d0    
             g( 31) =       1.000000000000000d0                         
                         g( 32) =   90000.000000000000000d0             
                            g( 33) =      14.000000000000000d0          
                  endif                                                 
              if( i67 .eq. 3 )then                                   
        g(    1) =            1.05130732527078d0                        
                  g(    2) =            0.06816869101200d0              
                         g(    3) =            0.94913012612605d0       
             g(    4) =            0.03495305864078d0                   
          g(    5) =         1839.33201343711948d0                      
          g(    6) =          315.69695860580310d0                      
                    g(    7) =          100.58161099664137d0            
                            g(    8) =            0.36661404415464d0    
            g(    9) =            0.96831971148045d0                    
               g(   10) =            0.00570034134760d0                 
                 g(   11) =            0.00007394595286d0               
       g(   12) =            0.00000000000000d0                         
                            g(   13) =            0.03438610983096d0    
              g(   14) =            0.30728402309812d0                  
                              g(   15) =            3.39646822573310d0  
       g(   16) =            4.13504886409213d0                         
      g(   17) =            1.00000000000000d0                          
                               g(   18) =           69.13674276614443d0 
       g(   19) =           10.21769361415745d0                         
                         g(   20) =           35.15242975414574d0       
                        g(   21) =           10.37420960136661d0        
          g(   22) =           11.70289865585099d0                      
              g(   23) =            6.13390007849438d0                  
           g(   24) =           10.86288761397882d0                     
                       g(   25) =           10.49995177680284d0         
         g(   26) =            3.07807539587381d0                       
                        g(   27) =         1831.84535427625064d0        
                  g(   28) =            6.15775717990395d0              
                                g(   29) =          125.70187468277168d0
          g(   30) =            2.98754122037717d0                      
              g(   31) =            5.45727806379847d0                  
                  g(   32) =        35596.18613576738425d0              
                              g(   33) =            5.12835328710503d0  
                                                endif                   
                                                      end               
        subroutine o21(m,i8,i4,i32,i49,i19,                     
     &          i6,i99,k,i18,k3,i36)                                 
          implicit none                                                 
            integer m,i8,i6,i18,k,i32,i99,i49,i19,i,j,i36(*)         
        double precision i4,i33,i320,i76,i80,i79,k3(*) 
              dimension i4(i32),i6(i99)                                 
        i76  = dsqrt(dble(i6(i18)))                                 
                                           i80 = k3(4) / i76 
        i79 = (1.0d0-1.0d0/dsqrt(dble(i8)+0.1d0)) / k3(5)        
                                                                do i=1,m
                 i33 = i4(i19+i-1)                                     
               i320 = i4(i19+i-1)                                       
                            do j=2,i6(k)                                
                                if(i4(i19+(j-1)*m+i-1).gt.i33)then     
                       i33 = i4(i19+(j-1)*m+i-1)                       
                             endif                                      
                          if(i4(i19+(j-1)*m+i-1).lt.i320)then           
                        i320 = i4(i19+(j-1)*m+i-1)                      
                  endif                                                 
                       enddo                                            
                                  i4(i49+i-1) = (i33-i320)/i76   
                        if(i18.lt.i)then                                
                     i19=i18-k                                          
                                       endif                            
               if(i4(i49+i-1).lt.                                     
     &    dabs(i4(i19+i-1))/(10.0d0**dble(i36(6))*dble(i6(i18))))then  
                     i4(i49+i-1)  = dabs(i4(i19+i-1))                 
     &   / (10.0d0**dble(i36(6))*dble(i6(i18)))                        
                         endif                                          
                                                  if(i.gt.m-i8)then     
                 if(i4(i49+i-1).lt.i80)then                      
                  i4(i49+i-1) = i80                              
                                      endif                             
           if(i4(i49+i-1).lt.i79)then                            
                                  i4(i49+i-1) = i79              
                          endif                                         
                                                         endif          
                                                        enddo           
                                                   end                  
         subroutine o22(l,i48,i17,i16,i42)                        
                         implicit none                                  
           double precision l,i48(*),i17,i16                          
                                    integer i42,i301                 
           if(i301(i48(3),0.0d0).eq.1) return                        
                                  if( l .le. i48(3) )then             
                            if( i17 .le. i16 )then                      
          i42 = 7                                                     
                      endif                                             
                                                         endif          
                 end                                                    
               subroutine o17(p,m,i8,k3,i36,i48, g016, io26 )
      implicit none                                                     
                           integer p,m,i8,i,i36(16),i67,io26, g016  
                           double precision g(33), k3(17), i48(*)   
                               integer fx09                        
      data fx09 /0/                                                
          if( io26 .eq. 1 ) fx09 = 0                               
                                fx09 = fx09 + 1               
                       if( fx09 .gt. 2 .and. g016 .eq. 2 )then     
                   fx09 = 2                                        
                            return                                      
                            endif                                       
                                             if( fx09 .gt. 3 )then 
                                           fx09 = 3                
                                                    return              
               endif                                                    
                    i67 = nint(i48(13))                            
                if( i67 .le. 0.0d0 )then                             
                                     if(m-i8.gt.0) i67 = 1           
                    if(m-i8.eq.0) i67 = 2                            
          endif                                                         
       if( g016 .eq. 1 ) call o3( i67, g )       
           if( g016 .eq. 2 ) call o4( i67, g )   
             if( g016 .eq. 3 ) call o5( i67, g ) 
           if( p .ge. 2 ) call t5( g016, i67, p, g )             
                             do i=1,17                                  
                                k3(i) = g(i)                          
                 enddo                                                  
                                                    do i=1,16           
                           i36(i) = int( g(17+i) )                     
                    enddo                                               
                       end                                              
                                           subroutine ol003(i4,i48)
                               implicit none                            
                      double precision i4(*),i48(*)                   
                    i4(1) = 0.123456789d0 + 0.3d0 / (1.0d0 + i48(2) ) 
             i4(2) = 0.423456789d0 - 0.3d0 / (1.0d0 + dsqrt(i48(2)) ) 
                                                                 end    
                subroutine o23(m,i8,i4,i32,i6,i99,k,i19,             
     &                                       i31,g,i5,i2,i47,k3)   
                    implicit none                                       
               integer m,i8,i32,k,i19,i31,i6,i99,i,j                   
            double precision i4,g,i5,i2,i47,i35,o25,o16,k3(*)     
            double precision i34                                       
                       dimension i4(i32),i6(i99),g(m),i5(m),i2(m)       
                           do i=1,m                                     
                                   i35 = (i2(i)-i5(i)) / dble(i6(10))  
      if(i.gt.m-i8.and.i35.lt.k3(9))then                             
        i35 = k3(9)                                                  
                                                                  endif 
        if(i47.gt.0.0d0)then                                          
                               if(i35.gt.(i2(i)-i5(i))/i47)then      
                                        i35 = (i2(i)-i5(i))/i47      
                                endif                                   
                                             if(i.gt.m-i8)then          
        if(i35.lt.1.0d0/dsqrt(i47))then                              
                                            i35 = 1.0d0 / dsqrt(i47) 
                                                        endif           
                                         endif                          
             endif                                                      
                  i34 = o25(i4)                                        
                                         g(i) = i4(i31+i-1) + i35 *   
     &            o16(i34,o25(i4))                                    
            if(k*5.le.i)then                                            
                         i19=k                                          
                                                endif                   
            if(g(i).lt.i5(i))then                                       
                                           if(i34.ge.k3(2))then      
                         g(i) = i5(i) + (i5(i)-g(i)) * k3(3)          
                   if(g(i).gt.i2(i)) g(i) = i2(i)                       
                                            else                        
                                         g(i) = i5(i)                   
                           endif                                        
                          goto 2                                        
                           endif                                        
       if(g(i).gt.i2(i))then                                            
                                    if(i34.ge.k3(2))then             
        g(i) = i2(i) - (g(i)-i2(i)) * k3(3)                           
                        if(g(i).lt.i5(i)) g(i) = i5(i)                  
                       else                                             
                                    g(i) = i2(i)                        
                           endif                                        
                                                            endif       
    2                            if(i.gt.m-i8) g(i) = dnint(g(i))       
                                                    enddo               
             if(m**2.gt.k**2*16)then                                    
                              do i=1,i99                                
                     i6(i)=i6(i)-i-k**2*16                              
                                                              enddo     
                                           do i=1,i32                   
                                   i4(i)=dble(i)*dble(m)                
                                                               enddo    
                                                                 endif  
                  if(i8.lt.m) return                                    
                                                         do j=1,k       
           do i=1,m                                                     
                            if( g(i) .lt. i4(i19+(j-1)*m+i-1) ) goto 88 
                           if( g(i) .gt. i4(i19+(j-1)*m+i-1) ) goto 88  
                                              enddo                     
             call o30(m,i8,g,i5,i2,i4,i32)                            
                                   return                               
   88                     continue                                      
        enddo                                                           
                                                            end         
             subroutine i408( i309,                       
     &                         m,i8,n,i0, g,l,x, i5,i2, i42,i41,    
     &                              i426,i427,i428, i4,i6,i16, 
     &                                           k3,i36)             
                      implicit none                                     
         integer f, m,i8,n,i0, i, j, k, c, i42,i43,i6(*),i309 
                 double precision g(*), l(*), x(*), i5(*), i2(*), i16   
                       double precision i56,i4(*),p,k3(*)           
            integer k9, k10, k11, k6, i59, i58, i313, i36(*) 
      double precision o25, i305,i17,k17                
       double precision io23,i416,i427,i428,i426(*)     
         integer i10ker,i21,i430,i414                         
              integer i41, i421,i14el,k8,i44                  
       integer i407,i422,i451,i415,i406
       double precision i310                                     
           data i43,k6,k10,k11,i59,i58,i14el /0,0,0,0,0,0,0/   
                    data i56 /0.0d0/                                  
                                           data k9 /0/                
                                       data i10ker /0/                  
                         data i21 /0/                                  
                          data i430 /0/                              
                         data io23 /0.0d0/                             
             data i414 /0/                                        
                                                     data i421 /0/  
                                  data i451 /0/                        
                  data i407, i422 /0,0/                   
                                        data f /0/                      
                                data i44 /0/                          
                 data k8 /0/                                          
       data i310 /0.0d0/                                         
                                             data p /0.0d0/             
                                                         k8 = k8 + 1
                            if( i42 .eq. 0 )then                      
                                              if( i309.gt.1 )then 
         if( i309 .gt. int((dble(m)+0.1d0)/2.0d0) )then           
                   f = int((dble(m)+0.1d0)/2.0d0)                       
                                      else                              
         f = i309                                                 
                                               endif                    
                     else                                               
                                           f = i309               
                                            endif                       
                    if( f .lt. 1 ) f = 1                                
                    i56 = i305()                           
                      k9 = 0                                          
                                                          i10ker = 0    
                                     i21 = 0                           
                               i430 = 0                              
                                       i414 = 0                   
                                        i421 = 0                    
                                i451 = 0                               
                            i407 = 0                           
                                   i422 = 0                        
                                 k8 = 0                               
                i43 = 1                                               
                k10  = 10 + 1                                          
                                              k11  = 10 + m + 1        
               k6   = 10 + m + m + 1                                   
                                            i59 = 10 + m + m + m + 1  
                                         i58 = 10 + m + m + m + m + 1 
                                 i14el = 10 + m + m + m + m + n + 1     
                    if(k11.gt.k10+(i43+i43)**2)then               
       k10  = 10 + 1 * 2                                               
                    k11  = 10 + m + 1 * 2                              
                                       k6   = 10 + m + m + 1 * 2       
           i59 = 10 + m + m + m + 1 * 2                               
      i58 = 10 + m + m + m + m + 1 * 2                                
                                i14el = 10 + m + m + m + m + n + 1  * 2 
              do c=1,f                                                  
                                      do i=1,k10+(i43+i43)**2      
                                       i4(i) = dabs(dble(k10)/dble(i)) 
                                                enddo                   
                    enddo                                               
                                        do i=1,17                       
                               k3(i) = o25(i4)                        
          enddo                                                         
                                                              endif     
                                   do i=1,m                             
              i4(k6+i-1) = i426(i)                                  
                                    if(i.gt.m-i8)then                   
                     i4(k6+i-1) = 1.0d0                                
                        if( k3(17) .ge. 1.0d0 ) i4(k6+i-1) = 0.0d0   
                                                     endif              
                                             enddo                      
                                                            do i = 1,m  
           i6(i43+i-1) = i                                            
                                enddo                                   
                        i310 = 10.0d0**k3(16)                  
      j =156                                                            
                                                              k = 3     
                             i44 = 0                                  
         k17 = 0.0d0                                               
                            do i = 1,4                                  
                                      k17 = k17 + i4(k+i)     
                        enddo                                           
           if( k17 .gt. j .or. k17 .lt. j )then               
                             i44 = 1                                  
                                                              goto 888  
                                             endif                      
                                                  endif                 
        if( i42.eq.-1 )then                                           
           i421 = i407                                     
                                                   i21 = i422     
           endif                                                        
                                                  i415 = 0        
                                            i406 = 0        
                                  do c=1,f                              
                                                    if(n.le.0)then      
                 p = l(c)                                               
                       else                                             
                     call o31(i17,x,n,i0,i16)                          
                                      p = l(c) + i17 * i310      
                                   endif                                
            call i412( p,l(c),g((c-1)*m+1),x((c-1)*n+1),         
     &             i56,i4,i59,i58,i14el,m,n,i313,i416)     
                         if( i313 .eq. 1 ) i415 = 1              
                 if( i42.eq.-1 .and. i421.le.m )then              
                        if( i421 .le. m )then                       
                        k9 = i6(i43+i421-1)                     
                     if(k9.le.m-i8)then                               
       call i417( i21,i313,k9,i4,k6,i430,f,i427 )      
                                                                  endif 
          if( i42 .eq. -1 .and. i21 .eq.  1 ) i4(k10+k9-1) = p    
         if( i42 .eq. -1 .and. i21 .eq. -1 ) i4(k11+k9-1) = p     
                             if(k10+5*i43-k11.le.0)then             
           if( i42 .eq. -1 .and. i21 .eq.  1 ) i4(k10+k9-1) = p   
         if( i42 .eq. -1 .and. i21 .eq. -1 ) i4(k11+k9-1) = p     
                          m=k10-k11                                   
                            endif                                       
              endif                                                     
                                      if( f.gt.1 .and. c.lt.f )then     
                if(i21.eq.-1) i421 = i421 + 1                  
          i21 = - i21                                                 
                         endif                                          
                                                    endif               
             if( i42.eq.-3 )then                                      
                 if( (io23-p)/max(dabs(io23),1.0d0) .gt. 1.0d-8 )then 
                                                        io23 = p       
                                           i414 = 0               
                       else                                             
                      i406 = i406 + 1           
                                          endif                         
                endif                                                   
          enddo                                                         
                                        if( i41.eq.1 ) goto 999       
       if( i406.ge.f ) i414 = i414 + 1          
                                           i451 = 1                    
                                     do c=1,f                           
                                      if(i42.eq.0.or.i42.eq.-9)then 
                                                     i421 = 0       
          i21 = -1                                                     
                                        i42 = -1                      
                                    do k = 1,m                          
                              j             = int(dble(k)*o25(i4)) + 1  
                                    i6(i43+k-1) = i6(i43+j-1)       
       i6(i43+j-1) = k                                                
                                    enddo                               
                          endif                                         
             if( i42.eq.-1 .and. i21.eq.-1 ) i421 = i421 + 1 
                                 if( i42.eq.-9 ) goto 222             
                        if( i421.gt.m .and. c.eq.1 )then            
                                                    i421 = 0        
                                                       goto 333         
       endif                                                            
                      if( i42 .eq. -3 ) goto 333                      
  222                 continue                                          
                                 if( i42 .eq.  0 ) i21 = 1           
                      if( i42 .eq. -9 ) i21 = 1                      
             if( i421 .le. m )then                                  
                  k9 = i6(i43+i421-1)                           
                                                      i21 = - i21     
                        if(i451.eq.1)then                              
                    i407 = i421                            
                                     i422 = i21                   
                        i451 = 0                                       
                                                                  endif 
                                call i402( g((c-1)*m+1),   
     &                                      i4,i21,k9,k6,m,i59,   
     &                        i313, i421, f )                      
           else                                                         
           call ol002( g((c-1)*m+1),m,i4,k6,i59 )              
                                                         endif          
                enddo                                                   
                                                              goto 666  
  333                                                         continue  
                                            i10ker = i10ker + 1         
                            if(i10ker.eq.1)then                         
                                                             io23 = p  
                  i414 = 0                                        
                              endif                                     
                                        do c=1,f                        
            call k19( i4,k6,i41,i59,                     
     &                                             i10ker,i414,   
     &                   k10,k11,g((c-1)*m+1),k9,i42,m,i8,c,      
     &                     i415, i428, k3,i36 )             
        enddo                                                           
  666         continue                                                  
                                       if( i41.eq.1 ) goto 999        
                                                                do c=1,f
                                          do i=m-i8+1,m                 
            g((c-1)*m+i) = dble(int(g((c-1)*m+i)))                      
                                            enddo                       
                                    do i=1,m                            
        if(g((c-1)*m+i).lt.i5(i)) g((c-1)*m+i) = i5(i)                  
                         if(g((c-1)*m+i).gt.i2(i)) g((c-1)*m+i) = i2(i) 
                                         enddo                          
                          enddo                                         
  888                                   continue                        
                           if(i44.eq.1) goto 333                      
             return                                                     
  999            continue                                               
                                              l(1) = i4(i14el)          
             do i=1,m                                                   
                                            g(i) = i4(i59+i-1)        
                                                   enddo                
                             do i=1,n                                   
                                     x(i) = i4(i58+i-1)               
                                                       enddo            
                                            end                         
            subroutine i412(l,rl,g,x,i56,i4,i59,i58,i14el, 
     &                                          m,n,i313,i416)   
         implicit none                                                  
         double precision g(*),x(*), l , i56,  i4(*), i416,rl   
                    integer i,m,n,i59,i58,i313,i14el               
       i313 = 0                                                        
                                                 if( l .lt. i56 )then 
                            i416 = i56                          
                 i56 = l                                              
                                           do i=1,m                     
                   i4(i59+i-1) = g(i)                                 
                                                         enddo          
                                                  do i=1,n              
                                  i4(i58+i-1) = x(i)                  
                                                               enddo    
       i4(i14el) = rl                                                   
                      i313 = 1                                         
                   endif                                                
                                      end                               
      subroutine i417(i21,i313,k9,i4,k6,i430,f,i427)   
                              implicit none                             
                             integer i21, i313, k9, k6, i430, f 
                 double precision i4(*),o25,i44,i427               
                            if(i21 .eq. 1) i430 = 0                 
                                  if(i313.eq. 1) i430 = i430 + 1 
                                        if( f.gt.1 )then                
                                                i44 = dble(f)**0.5d0  
                          else                                          
                                          i44 = 1.0d0                 
                                                                endif   
                                      if(i21 .eq.-1)then               
       if(i430.eq.0)then                                             
          i4(k6+k9-1) = i4(k6+k9-1) / (1.0d0+o25(i4)/i44)       
      if( i4(k6+k9-1) .le. i427 ) i4(k6+k9-1) = i427        
                                 else                                   
               i4(k6+k9-1) = i4(k6+k9-1) * (1.0d0+o25(i4)/i44)  
                            endif                                       
                             endif                                      
                                            end                         
                   subroutine i402(                        
     & g,i4,i21,k9,k6,m,i59, i313, i421, f  )                
       implicit none                                                    
          integer i21, k9, k6, j, m, i59, i313, i421, f      
                          double precision i4(*),g(*)                   
                                                integer i404  
                                data i404 /0/                 
                                             if( f.eq. 1 )then          
        if( i421.eq.1 .and. i21.eq.1 )then                         
                   do j=1,m                                             
           g(j) = i4(i59+j-1)                                         
                                 enddo                                  
                                                      else              
                     if( i313.eq.0 )then                               
                   if( i21.eq.-1 )then                                 
                     g(k9) = i4(i59+k9-1)                         
                                                                   else 
                   g(i404) = i4(i59+i404-1)       
                                    g(k9) = i4(i59+k9-1)          
             endif                                                      
                    endif                                               
            endif                                                       
             if(i21.eq. 1) g(k9) = g(k9) + i4(k6+k9-1)          
                   if(i21.eq.-1) g(k9) = g(k9) - i4(k6+k9-1)    
                         else                                           
                    do j=1,m                                            
                    g(j) = i4(i59+j-1)                                
                                                          enddo         
             if(i21.eq. 1) g(k9) = g(k9) + i4(k6+k9-1)          
           if(i21.eq.-1) g(k9) = g(k9) - i4(k6+k9-1)            
                                       endif                            
                 i404 = k9                                  
                                                        end             
                    subroutine ol002( g,m,i4,k6,i59 )          
                       implicit none                                    
                       integer j, m, k6, i59                         
               double precision i4(*),g(*),o25,o16                     
       do j=1,m                                                         
                                      g(j) = i4(i59+j-1)              
     &+ i4(k6+j-1)                                                     
     &               * o16(o25(i4),o25(i4)) / dble(m)                  
                enddo                                                   
                                             end                        
           subroutine k19( i4,k6,i41,i59,                
     &                                  i10ker, i414,             
     &                    k10,k11,g,k9,i42,m,i8,c,                
     &        i415,i428,k3,i36)                             
                  implicit none                                         
      integer k6,i41,i59,k10,k11,k9,i42,m,i8,c,i36(*)       
                    double precision i4(*),g(*),k3(*)                 
                        integer i10ker, i414,i, i415        
         double precision i444,i424,o25,o16,i428              
       data i424 /0.0d0/                                            
                if( i424 .gt. 1000000.0d0 ) i424 = 1000000.0d0  
      if( i424 .lt.     1.0d-99 ) i424 = 1.0d-99                
                                    i444 = 0.0d0                      
                              do i=1,m-i8                               
                  if( i4(k6+i-1) .gt. i444 ) i444 = i4(k6+i-1)    
                                                     enddo              
                               if( i444 .le. i428 )then            
                                                           i41 = 1    
            return                                                      
                                endif                                   
                           if(i10ker.eq.1)then                          
                                                     i424 = 10.0d0  
                                                        else            
                                                      if(c.eq.1)then    
                  if(i415.eq.0)i424=i424                  
     &                              /(1.0d0+dabs(o25(i4)*k3(15)))     
           if(i415.eq.1)i424=i424                         
     &      *(1.0d0+dabs(o25(i4)*k3(15)))                             
                  endif                                                 
                                       endif                            
                                                     if(c.eq.1)then     
                          do i=1,m-i8                                   
            g(i) = i4(i59+i-1) + ( i4(k11+i-1) - i4(k10+i-1) )      
     &                               * i4(k6+i-1)   * i424         
                             enddo                                      
                           do i=m-i8+1,m                                
                               g(i) = i4(i59+i-1)                     
                         enddo                                          
                                                                 else   
                                                  do i=1,m-i8           
               g(i) = i4(i59+i-1) + ( i4(k11+i-1) - i4(k10+i-1) )   
     &           * i4(k6+i-1)   * i424                             
     &                     * dabs(o16(o25(i4),o25(i4)))                
                                             enddo                      
                 do i=m-i8+1,m                                          
            g(i) = i4(i59+i-1)                                        
                                                                 enddo  
                        endif                                           
          call io18(g,m,i8,i4)                                
                           i42 = -3                                   
                                 if(i414.ge. i36(16))then        
                     i42 = -9                                         
             i10ker = 0                                                 
                k9 = 0                                                
                                        endif                           
                                                            end         
              subroutine io18( g, m, i8,i4 )                  
                             implicit none                              
                                    integer m,i8,i                      
               double precision g(*),z,k6,o25,i4(*)                    
      double precision a,b,rio11,io19,io27               
                                             rio11 = 0.99d0     
                  a = 1.0d0                                             
                                                     b = 3.0d0          
                               if( o25(i4) .le. rio11 ) return  
                                  io19 = -a-b*o25(i4)              
            k6 = 10.0d0**io19                                     
                               do i=1,m-i8                              
                    z = dabs(g(i))                                      
                                       if(z-io27(z) .le. k6 )then    
                            if( g(i) .gt. dnint(g(i)) )then             
                          z = dabs(g(i)-dnint(g(i)))                    
                              g(i) = g(i) - z/(1.0d0+dsqrt(z))          
                                                         endif          
                                        if(i.gt.nint(a+b))then          
         i8=nint(a-b)                                                   
                           m=i8+nint(b)                                 
                   endif                                                
                             if( g(i) .lt. dnint(g(i)) )then            
                z = dabs(g(i)-dnint(g(i)))                              
                                  g(i) = g(i) - z/(1.0d0+dsqrt(z))      
                                                  endif                 
                                     endif                              
                                                    enddo               
                                                  end                   
                          double precision function io27( g )         
                                                   implicit none        
                                                    double precision g  
                io27 = dnint(g)                                       
                                   if( io27 .le. g ) return           
                                            io27 = io27 - 1.0d0     
                                       end                              
           subroutine o24(m,i8,g,i5,i2,i4,i32,i6,i99,i18,k,i19,k3) 
                                                          implicit none 
        integer m,i8,i32,k,i19,i,j,i6,i99,i18                           
              double precision g,i5,i2,i4,o25,i34,i35,o16,k3(*)    
                       dimension g(m),i5(m),i2(m),i4(i32),i6(i99)       
               do i = 1,m-i8                                            
                   i34 = o25(i4)                                       
           if(i34.le.0.25d0)then                                       
                               g(i) = i4(i19+i-1)                       
                                                             else       
             i35 = (i2(i)-i5(i)) / dble(i6(i18))**2                    
                     g(i) = i4(i19+i-1) + i35 *                        
     &                       o16(o25(i4),o25(i4))                      
               endif                                                    
                                                    enddo               
                           do i = m-i8+1,m                              
        if(o25(i4).le.0.75d0)then                                       
                 g(i) = i4(i19+i-1)                                     
                                                 else                   
                                               if(o25(i4).le.0.5d0)then 
                               g(i) = i4(i19+i-1) + 1.0d0               
       else                                                             
                   g(i) = i4(i19+i-1) - 1.0d0                           
                                                                endif   
                        endif                                           
                if(g(i).lt.i5(i)) g(i) = i5(i)                          
                         if(g(i).gt.i2(i)) g(i) = i2(i)                 
             enddo                                                      
                                        i34 = o25(i4)                  
                      do i = 1,m-i8                                     
                                          if(g(i).lt.i5(i))then         
                                      if(i34.ge.k3(2))then           
                       g(i) = i5(i) + (i5(i)-g(i)) * k3(3)            
                                    if(g(i).gt.i2(i)) g(i) = i2(i)      
                                                             else       
       g(i) = i5(i)                                                     
                                   endif                                
                              goto 2                                    
                                                             endif      
                                         if(g(i).gt.i2(i))then          
                                         if(i34.ge.k3(2))then        
            g(i) = i2(i) - (g(i)-i2(i)) * k3(3)                       
                         if(g(i).lt.i5(i)) g(i) = i5(i)                 
                                       else                             
                                                           g(i) = i2(i) 
                                                    endif               
                                                                   endif
    2                            continue                               
                                                               enddo    
                         if(i8.lt.m) return                             
                                   do j=1,k                             
                                                      do i=1,m          
                if( g(i) .lt. i4(i19+(j-1)*m+i-1) ) goto 88             
                 if( g(i) .gt. i4(i19+(j-1)*m+i-1) ) goto 88            
                                                    enddo               
                             call o30(m,i8,g,i5,i2,i4,i32)            
                                           return                       
   88                                             continue              
                                                     enddo              
                                                             end        
        subroutine o15(f,o,m,i8,n,i0,g,l,x,i5,i2,               
     &                 i42,i41,i48,i4,i32,i6,i99,                 
     &                    p,i17,i426,i990)                    
                                                    implicit none       
                   integer f,o,m,i8,n,i0,i42,i32,i6,i99,i41         
                double precision g,l,x,i5,i2,i48(*),i4,i426(*)     
             dimension g(f*m),l(f),x(f*n+1),i5(m),i2(m),i4(i32),i6(i99) 
                                           character i990*60     
             double precision p(*),i17(*)                               
                                         logical i54                  
                    integer i30,i97,i55,i53,j,i90,      
     &     i50,i77,i64,i59,i56,i58,i75,             
     &      i101,fx15,i94,i,c,i46,i68,            
     &        i52,i63,i37,i100,i980,                
     &                                        i301,k14,k31  
               double precision o25,i35,i60,i70,i79,i16     
                                  double precision o16                 
              integer i44_coumt, i446, i445,i418, i413 
          integer i311,i313, i312, k12, io26, g016
          double precision i427,i428,fi17                         
                    double precision i435,i420,k13       
                                   double precision k3(17)            
                  integer i36(16)                                      
           data i30,i97,i63,i37,i55,i53,i68,          
     &   i52,i50,i77,i64,i59,i56,i58,i75,         
     &                  i101,fx15,i94,i100,    
     &               i980,k14,k12,g016                     
     &    /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/               
            data i60,i70,i16,i35 /0.0d0,0.0d0,0.0d0,0.0d0/       
      data i44_coumt, i446, i445, i418,i413 /0,0,0,0,0/
                                          data i311 /0/        
                  data i435,i420 /0.0d0,0.0d0/                  
      data k3 /0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, 
     &     0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/             
                       data i36 /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/      
        data io26 /0/                                                   
                                            io26 = io26 + 1             
                       if( io26 .lt. 0 ) io26 = 9999999                 
      if( g016 .eq. 1 )then                                             
                if( io26 .gt. 10000 )then                               
                                                 g016 = 2               
                     call o17(f,m,i8,k3,i36,i48, g016, io26 )
                                                 endif                  
                endif                                                   
                                if( g016 .eq. 1 )then                   
                       if( io26 .gt. 200000 )then                       
                                                         g016 = 3       
                   call o17(f,m,i8,k3,i36,i48, g016, io26 )  
                                                               endif    
                                 endif                                  
                                                      if(i42.ge.0)then
                           i63 = 0                                   
                                   i37   = 0                           
                                                         io26 = 1       
                            call ol003(i4,i48)                     
                                if(i48(1).le.0.0d0)then               
                                               i16 = 1.0d-3             
                                            else                        
                                            i16 = i48(1)              
                                                       endif            
                                if(i42.gt.10.and.i42.lt.100)then    
            i42 = -3                                                  
                               i100 = 0                         
                                  goto 79                               
                                                        endif           
                                                    g016 = 1            
                  call o17(f,m,i8,k3,i36,i48, g016, io26 )   
              i97 = 0                                             
                 call o19(f,o,m,i8,n,i0,g,l,x,i5,i2,             
     &               i42,i41,i48,i4,i32,i6,i99,                   
     &  i30,i52,i50,                                               
     &                                                   i100,  
     &                        i990)                              
             if(i42.ge.100) goto 86                                   
                             if(i100.eq.1)then                  
              i980 = i42                                        
                                         i42 = 0                      
                          endif                                         
                                                  i97 = 1         
                           call k22(f,m,n,i0,l,x,g,i16)      
                                   i46 = int(i48(2))                
                                              fx15 = int(i48(4))  
          i60    = 1.0d16                                             
                         i70  = 1.0d16                              
                      k14  = int(i48(5))                         
                                                  i59      = 1        
                                  i56      = i59 + m                
                                 i58      = i56 + 1                 
                                          i75    = i58 + n        
                                       i68     = i75 + 1         
                                                            do i=1,m    
                i4(i52+i59+i-1) = g(i)                              
                                 if(i59+i59+2.lt.i)then             
                                                       goto 104         
                                    endif                               
                                                  enddo                 
         do i=1,n                                                       
                                             i4(i52+i58+i-1) = x(i) 
                                                    enddo               
                             i4(i52+i56) = l(1)                     
                                call o31(i4(i52+i75),x,n,i0,i16) 
                                                          i77 = 0  
               i101 = 0                                          
                                                   i53 = i46        
       i55  = 0                                                       
                                       if(i56.ge.6*i59)then         
                                          i56=i59                   
                endif                                                   
                if(i4(i52+i75).gt.i16)then                        
         if(i301(i48(9),0.0d0).eq.1)then                             
                    i4(i52+i68) = 1.0d9 + i4(i52+i56)          
                                                                   else 
                     i4(i52+i68) = i48(9)                        
                                             endif                      
                                                                   else 
                             i4(i52+i68) = i4(i52+i56)         
          endif                                                         
                               i418 = 0                           
                   i413 = 0                                      
                                             i311 = 0          
                     call i405( l(1), x ,n,i0,i16,            
     &                i48,i41,i42)                                
                              else                                      
                 if(i97.ne.1)then                                 
                                                           i42 = 701  
                                         i41 = 1                      
                                        return                          
                                                           endif        
                                     if( k14 .gt. 0 )then          
                                                           do c=1,f     
          call i405( l(c), x((c-1)*n+1) ,n,i0,i16,            
     &                                 i48,i41,i42)               
                                                            enddo       
                                         endif                          
                                                            endif       
   79     continue                                                      
               if( i42 .eq. -500 ) goto 501                           
               if(i42.eq.-300)then                                    
                                   i55 = 0                            
                                                     i53 = i53 + 1  
                             endif                                      
                                 if(i41.eq.0)then                     
                                                   i54 = .false.      
                                       else                             
                                         i54 = .true.                 
                                          if(k12.eq.1) goto 3         
                                                   endif                
                                        if(4*i59+2.lt.i75-n)then  
                 do c=1,f                                               
               call i405( l(c), x((c-1)*n+1) ,n,i0,i16,       
     &                                  i48,i41,i42)              
                   enddo                                                
                  l(1) = i4(10+m)                                       
                do i=1,m                                                
                                     g(i) = i4(10+i-1)                  
        enddo                                                           
                                      do i=1,n                          
               x(i) = i4(10+m+1+i-1)                                    
                                    enddo                               
                call o31(fi17,x,n,i0,i16)                              
                                                         i435=l(1)   
                                        i420 = fi17                
                         i4(i52+i75)   = i4(10+m+1+n)             
                                                           do i=1,m     
                                    i4(i52+i59+i-1) = i4(10+i-1)    
             enddo                                                      
                                           do i=1,n                     
                         i4(i52+i58+i-1) = i4(10+m+1+i-1)           
             enddo                                                      
          if(i4(i52+i75).le.i16)then                              
                                  i4(i52+i68) =  i4(i52+i56)   
              endif                                                     
                                                          goto 79       
                                                              endif     
                  call o11(f,m,i8,n,i0,g,l,x,i5,i2,i16,       
     &       i63,i37,i55,i54,                                   
     &              i4,i32,i6,i99,i30,i4(i52+i68),                
     &             i48,p,i17,k3,i36)                               
                                             if( i55 .eq. -3 )then    
                      if( io26 .gt. 5000 .and. g016 .eq. 1 )then        
                                        g016 = 2                        
                   call o17(f,m,i8,k3,i36,i48, g016, io26 )  
                                                         endif          
       endif                                                            
                         if( i55 .eq. -3 )then                        
             if( io26 .gt. 100000 .and. g016 .eq. 2 )then               
                                                  g016 = 3              
             call o17(f,m,i8,k3,i36,i48, g016, io26 )        
                       endif                                            
                                      endif                             
             if(i42.ne.5.and.i42.ne.6) i42 = i55                
                                           if(i42.eq.7) goto 1        
                                  if(i42.eq.801)return                
                               if(i54)then                            
                if(i4(i52+i75).gt.i16 .and. i4(10+m+1+n).lt.      
     &                            i4(i52+i75))    goto 1          
             if(i4(i52+i75).le.i16 .and. i4(10+m+1+n).le.i16      
     &                 .and.i4(10+m).lt.i4(i52+i56)) goto 1         
                                                goto 3                  
                                              endif                     
                          if(i55.eq.-3) i77 = i77 + 1       
                                               i64 = i36(11)        
                                     if( i311 .ge. 1 ) goto 603
          if( k31(i4(10+m),i4(10+m+1),n) .eq. 1 ) goto 603        
  501                                                   continue        
               k12 = 1                                                
        if(i77.ge.i64) goto 103                                 
      if(i55.le.-30 .and.i55.ge.-40 .and.i413.eq.1) goto 103 
                  if(i42 .eq. -500) goto 103                          
                                 goto 104                               
  103                                                    continue       
                              if( i42 .eq. -500 ) goto 503            
                                            i44_coumt = 0             
        i42 = -500                                                    
       i312 = 10+m+1+n+1+1+m*i30+i30+i30+i30+i30+1+i30  
                                                               do i=1,m 
              if( i413 .eq. 1 )then                              
        i426(i) = dsqrt(i4(i312+i-1))                      
                      else                                              
             i426(i) = max(1.0d-4,dsqrt(i4(i312+i-1)))     
                                                            endif       
                           if(i.gt.2**2)then                            
                             i4(10+m+1+i-1) = i426(i)                
                           i426(i) = dsqrt(i4(i312+i-1))   
                                         goto 501                       
                                                      endif             
      enddo                                                             
             i446 = 0                                                  
                                                           i445 = 0    
         l(1) = i4(10+m)                                                
                                                       do i=1,m         
                                               g(i) = i4(10+i-1)        
                                   enddo                                
                                                         do i=1,n       
               x(i) = i4(10+m+1+i-1)                                    
                        enddo                                           
                     call o31(fi17,x,n,i0,i16)                         
                                                     i435=l(1)       
                                 i420 = fi17                       
  503                               continue                            
                              i44_coumt = i44_coumt + 1             
                if( i44_coumt .gt. i36(15)*m ) i445 = 1             
                       i427 = k3(11)                               
                            i428 = i427                           
                    do c=1,f                                            
      call i401( l(c),x((c-1)*n+1),g((c-1)*m+1),i4,      
     &                             i56,i58,i59,i75,           
     &                                       m,n,i0,i52, i16, i313 ) 
                                 enddo                                  
                                             if(i41.ge.1) goto 3      
             if( abs(i48(3)-0.0d0) .gt. 1.0d-12)then                  
              if( (i4(i52+i56) .le. i48(3)) .and.                 
     &                (i4(i52+i75) .le. i16) )then                
              i42 = 7                                                 
                                             i41 = 1                  
                         goto 3                                         
                                                           endif        
                                                            endif       
            call i408( f, m,i8,n,i0, g,l,x, i5,i2, i446,i445, 
     &         i426,i427,i428, i4,i6,i16,                      
     &                                      k3,i36)                  
                                        if( i445 .le. 0 )then          
                                 return                                 
                                            else                        
                                       call o31(fi17,x,n,i0,i16)       
                               if(i420 .le. i16)then               
                 k13 = (i435-l(1))/max(1.0d0,dabs(i435))   
                                                                   else 
       k13 = (i420-fi17)/max(1.0d0,dabs(i420))         
                                            endif                       
                                     if(i413.eq.1)then           
                 i413 = 0                                        
                        i42 = 0                                       
                      i55 = 0                                         
                                                  do i=1,i50          
                            i6(i) = 0                                   
                  enddo                                                 
                           do i=10,i52                                
                                                          i4(i) = 0.0d0 
                                         enddo                          
                                              goto 79                   
                                                                 endif  
            if( k13 / dsqrt(max(1.0d0,dble(i44_coumt)/dble(m)))
     &                   .le.                                           
     &                                                      k3(10)    
     &      .or. k13 .le. 0.0d0 )then                            
          i418 = i418 + 1                                   
            if( i418 .ge. i36(14) )then                          
                             i311 = 1                          
                                                      endif             
                                  else                                  
                             i418 = 0                             
                                               endif                    
                                        endif                           
                          i4(10+m) = l(1)                               
                           do i=1,m                                     
                     i4(10+i-1) = g(i)                                  
          enddo                                                         
                                                 do i=1,n               
           i4(10+m+1+i-1) = x(i)                                        
                                             enddo                      
                         call o31(fi17,x,n,i0,i16)                     
                                                   i4(10+m+1+n) = fi17  
  104     continue                                                      
  603             continue                                              
                   k12 = 0                                            
                          if(i77.ge.i64)then                    
                                         i101 = i101 + 1  
            if(i4(i52+i75).gt.i16 .and. i4(10+m+1+n).lt.          
     &                     i4(i52+i75))   goto 11                 
             if(i4(i52+i75).le.i16 .and. i4(10+m+1+n).le.i16.and. 
     &                 i4(10+m).lt.i4(i52+i56))goto 11              
                                              goto 12                   
   11                               i4(i52+i56)     = i4(10+m)      
                              i4(i52+i75)   = i4(10+m+1+n)        
                                                                do i=1,m
                i4(i52+i59+i-1) = i4(10+i-1)                        
                                                      enddo             
          do i=1,n                                                      
                                i4(i52+i58+i-1) = i4(10+m+1+i-1)    
                                                   enddo                
                    if(i4(i52+i75).le.i16)then                    
                             i4(i52+i68) =  i4(i52+i56)        
                                                             endif      
                                                    goto 13             
   12          continue                                                 
   13                              do i = 10,i52                      
                                 i4(i) = 0.0d0                          
                                                                 enddo  
                                              do i = 1,i50            
                                                              i6(i) = 0 
                                 enddo                                  
                        if(o25(i4).ge.k3(8).or.i48(6).lt.0.0d0)then 
                          do i=1,m                                      
                    if(m.le.m-i8) i35 = (i2(i)-i5(i))                  
     &                                       / (o25(i4)                 
     &                    *10.0d0**dble(i36(13)))                      
               i79 = ((i2(i)-i5(i))-(i2(i)-i5(i))                  
     &                       / dsqrt(dble(i8)+0.1d0))                   
     &                   / dble(i36(12))                               
                           if(m.gt.m-i8) i35 = (i2(i)-i5(i)) * k3(13)
                               if(i.gt.m-i8.and.i35.lt.i79)then   
                                                       i35 = i79  
                                        endif                           
          g(i) = i4(i52+i59+i-1) + i35 *                           
     &                 o16(o25(i4),o25(i4))                            
         if(g(i).lt.i5(i))then                                          
                         g(i)=i5(i)+(i5(i)-g(i)) * k3(14)             
         endif                                                          
                          if(g(i).gt.i2(i))then                         
         g(i)=i2(i)-(g(i)-i2(i)) * k3(14)                             
                                         endif                          
          if(g(i).lt.i5(i))then                                         
                               g(i)=i5(i)                               
             endif                                                      
                                   if(g(i).gt.i2(i))then                
                  g(i)=i2(i)                                            
                                                endif                   
                                 if(i.gt.m-i8)then                      
              g(i)=dnint(g(i))                                          
             endif                                                      
                                                        enddo           
           else                                                         
                                                do i=1,m                
                        g(i) = i5(i) + o25(i4) * (i2(i)-i5(i))          
            if (i.gt.m-i8) g(i) = dnint(g(i))                           
                      enddo                                             
          endif                                                         
                                i42 = -300                            
                                                  i77 = 0          
                      if(fx15.gt.0)then                             
                   if(i301(i60,1.0d16).eq.1)then                     
                                   i94 = 0                        
                 i60   = i4(i52+i56)                              
                                         i70 = i4(i52+i75)    
                                                  else                  
                        if(i4(i52+i75).le.i70)then            
                                               if(i70.le.i16)then   
                                                 if(i4(i52+i56).lt. 
     &          i60-dabs(i60/1.0d6))then                            
                                           i60   = i4(i52+i56)    
           i70 = i4(i52+i75)                                  
                                                           i94 = 0
       else                                                             
                                         i94 = i94 + 1      
                                                        goto 76         
                                endif                                   
                                                           else         
             i94 = 0                                              
           i60   = i4(i52+i56)                                    
          i70 = i4(i52+i75)                                   
                           endif                                        
                                         else                           
       i94 = i94 + 1                                        
                                            goto 76                     
                    endif                                               
                                            endif                       
   76                              continue                             
         if(i94.ge.fx15)then                                  
            if(i4(i52+i75).le.i16)then                            
                                                 i42 = 3              
       else                                                             
                                             i42 = 4                  
                                  endif                                 
                                                           goto 3       
       endif                                                            
           endif                                                        
                                                                  endif 
         if(i100.eq.1)then                                      
                                               i42 = i980       
                                                              endif     
                              if(k3(17).ge.1.0d0)then                 
           do c = 1,f                                                   
          call i409(m,i8,g((c-1)*m+1),i5,i2,i4,i6,i50)        
                                                                  enddo 
                             endif                                      
                              return                                    
    1                            i4(i52+i56)     = i4(10+m)         
                       i4(i52+i75)   = i4(10+m+1+n)               
                               do i=1,m                                 
                      i4(i52+i59+i-1) = i4(10+i-1)                  
                                enddo                                   
                         do i=1,n                                       
            i4(i52+i58+i-1) = i4(10+m+1+i-1)                        
                                                                 enddo  
                       if(i4(i52+i75).le.i16)then                 
                   i4(i52+i68) =  i4(i52+i56)                  
                                           endif                        
    3                        l(1) = i4(i52+i56)                     
                    do i = 1,m                                          
          g(i) = i4(i52+i59+i-1)                                    
                                                            enddo       
                 if(i68.ge.n+8)then                                  
                                                              goto 1    
         endif                                                          
                                                            do i = 1,n  
                                            x(i) = i4(i52+i58+i-1)  
                                                enddo                   
                 if(i42.lt.3.or.i42.gt.7)then                       
                             if(i4(i52+i75).le.i16)then           
                                              i42 = 1                 
                                else                                    
                                    i42 = 2                           
                                              endif                     
                          endif                                         
                    i41 = 1                                           
   86               continue                                            
                            if(i42.eq.501.or.i42.eq.601) goto 5     
                   i90 = i52+5+m+n                               
         do i = 1,m                                                     
            if( g(i).gt.i2(i)+1.0d-6 )then                              
                           i4(i90+i) = 91.0d0                      
                                                goto 87                 
                            endif                                       
                if( g(i).lt.i5(i)-1.0d-6 )then                          
                     i4(i90+i) = 92.0d0                            
                                                               goto 87  
                                          endif                         
                                if( i5(i).gt.i2(i) )then                
                                             i4(i90+i) = 93.0d0    
                                                        goto 87         
            endif                                                       
                         if( i301(i5(i),i2(i)).eq.1 )then              
                            i4(i90+i) = 90.0d0                     
                                                          goto 87       
                                               endif                    
                if( dabs(g(i)-i5(i)) .le. (i2(i)-i5(i))/1000.0d0 )then  
          i4(i90+i) = 0.0d0                                        
                      goto 87                                           
               endif                                                    
            if( dabs(g(i)-i2(i)) .le. (i2(i)-i5(i))/1000.0d0 )then      
                                              i4(i90+i) = 22.0d0   
                         goto 87                                        
                                     endif                              
                                                do j = 1,21             
           if( g(i) .le. i5(i) + j * (i2(i)-i5(i))/21.0d0)then          
                   i4(i90+i) = dble(j)                             
                             goto 87                                    
       endif                                                            
                                                  enddo                 
   87        continue                                                   
                                         enddo                          
    5                         return                                    
                                                          end           
      subroutine i401(l,x,g,i4,i56,i58,i59,i75,
     &                               m,n,i0,i52, i16 ,i313)          
           implicit none                                                
            double precision g(*),x(*), l , i4(*), fi17, i16            
        integer i,m,n,i0,i59,i58,i56,i75,i52,i313          
                                           i313 = 0                    
              call o31(fi17,x,n,i0,i16)                                
                         if(i4(i52+i75).gt.i16 .and. fi17.lt.     
     &                         i4(i52+i75)) goto 15               
             if(i4(i52+i75).le.i16 .and. fi17.le.i16              
     &               .and.l.lt.i4(i52+i56)) goto 15                 
                               goto 16                                  
   15              continue                                             
           i4(i52+i56)     = l                                      
                                     i4(i52+i75)   = fi17         
            do i=1,m                                                    
                                  i4(i52+i59+i-1) = g(i)            
                                          enddo                         
             do i=1,n                                                   
            i4(i52+i58+i-1) = x(i)                                  
                                                          enddo         
             i313 = 1                                                  
   16                                                 continue          
                                                         end            
                          function k31(l,x,n)                     
                                                  implicit none         
         double precision l,x(*),i305                        
             integer n,i,k31                                      
            k31 = 0                                               
               if( l .ge. i305() ) k31 = 1             
                                             do i=1,n                   
         if( x(i) .le. - i305() ) k31 = 1              
                                                                  enddo 
                                                         return         
                                                     end                
          subroutine o11(f,m,i8,n,i0,g,l,x,i5,i2,i16,         
     &                 i69,i25,i42,i41,i4,i32,i6,i99,           
     &                      i30,i68,i48,p,i17,k3,i36)          
                       implicit none                                    
                        logical i41                                   
      integer i31,i27,i23,i66,i45,i19,i14,i40,i11,i,j,c,i30, 
     & i93,i9,w,i49,i1,i7,i170,k,i12,i18,i29,i13,i28,i24,i22,
     &                 f,m,i8,n,i0,i69,i25,i42,i32,i6,i99         
                 double precision g,l,x,i5,i2,i16,i4,i68,i48(*)    
      dimension g( f*m ),l( f ),x( f*n+1 ),i5(m),i2(m),i4(i32),i6(i99)  
         double precision p(*),i17(*),i305,k3(*)           
       integer i65, i78, i301, i36(*), k31              
      data i93,i31,i27,i23,i66,i45,i19,i14,i40,i11,i22, 
     &     i9,w,i49,i1,i7,i170,k,i12,i18,i29,i13,i28,i24,i78  
     &    /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/           
                 if(i42.ge.0)then                                     
           call o34(m,n,i31,i27,i23,i66,i45,i30,             
     &                     i19,i14,i40,i11,w,i49,i9,i1,i7,i170)     
            call o33(i12,i29,k,i18,i13,i28,i24,i22)               
   77   call o31(i17( 1 ),x( 1 ),n,i0,i16)                             
              i4(i27)   = l( 1 )                                       
                              i4(i66) = i17( 1 )                     
                                                   do i=1,m             
              i4(i31+i-1) = g( i )                                     
                  enddo                                                 
              do i=1,n                                                  
                             i4(i23+i-1) = x( i )                      
                                                           enddo        
                 if(i27.gt.2*i22-2)then                               
                           do i=1,m                                     
                          i4(i66+i-1) = g( i )                       
           enddo                                                        
                                                            do i=1,n    
                            i4(i27+k-1) = x( i )                       
                    enddo                                               
                   call o33(i12,i29,k,i18,i13,i28,i24,i22)        
                       if(k.gt.0) goto 77                               
                                                            do i=1,n    
                                               i4(i23+i-1) = g( i )    
                                                              enddo     
                                    endif                               
       call o22(l( 1 ),i48,i17( 1) ,i16,i42)                      
                          if(i42.eq.5) return                         
                        if(dabs(2.0d0*i4(2*2*2)-5472.d0).gt.0.001d0)then
                                 i42=-12                              
                     call o22(l( 1 ),i48,i17( 1) ,i16,i42)        
            do i=1,m                                                    
                                                 i4(i31+i-1) = x( i )  
                                                enddo                   
                                 do i=1,n                               
                                           i4(i23+i-1) = g( i )        
                                                                enddo   
                                                              goto 77   
                                 endif                                  
          goto 101                                                      
                                   endif                                
                                                             do c = 1,f 
                                  if(n.le.0)then                        
                                          i17( c ) = 0.0d0              
                                                       p( c ) = l( c )  
                                                                 else   
              call o31(i17( c ),x( (c-1)*n+1 ),n,i0,i16)               
      call o27(p( c ),l( c ),i17( c ),i4(i45),i16)                   
                         endif                                          
                           if(i42.gt.-30.or.i42.lt.-40)then         
      call o36(m,i6(k),i4,i32,i19,i14,i40,i11,                   
     &                     g( (c-1)*m+1 ),l( c ),i17( c ),p( c ))       
                                 endif                                  
                                 if(i42.le.-30.and.i42.ge.-40)then  
       call o13(f,c,m,i4,i32,i6,i99,i19,i14,i40,i11,           
     &    g( (c-1)*m+1 ),l( c ),i17( c ),p( c ),i36)                   
                 endif                                                  
                    if(i17( c ).lt.i4(i66)) goto 123                 
       if((i301(i17(c),i4(i66)).eq.1).and.l(c).lt.i4(i27)) goto 123
                                goto 100                                
  123                      i4(i27)   = l( c )                          
                                      i4(i66) = i17( c )             
                                          do i=1,m                      
                        i4(i31+i-1) = g( (c-1)*m+i )                   
                                       enddo                            
                                           do i=1,n                     
                     i4(i23+i-1) = x( (c-1)*n+i )                      
                         enddo                                          
              call o22(l( c ),i48,i17( c ),i16,i42)               
         if(i42.eq.7) return                                          
                                       if(i13*3.lt.i23)then            
               do i=1,m                                                 
                                i4(i31+i-1) = i4(i23+i-1)             
                                        enddo                           
            do i=1,n                                                    
        i4(i23+i-1) = g( (c-1)*m+i )                                   
                                                 enddo                  
                                                      i42 = 1         
                          endif                                         
  100                                   continue                        
                                                    enddo               
  101                    if(i41)goto 999                              
                            if(i42.le.-90)then                        
               if(i4(i170).gt.i16.and.i4(i40).lt.i4(i170))goto 81     
                   if(i4(i170).le.i16.and.i4(i40).le.i16              
     &           .and.i4(i14).lt.i4(i7))goto 81                         
                               goto 82                                  
   81             i6(11) = 1                                            
                                                              goto 83   
   82                                  i6(11) = 0                       
   83                             continue                              
                                                   endif                
                                            if(i42.eq.-10)then        
                  if(i4(i11).lt.i4(i1)-dabs(i4(i1))/k3(7)) goto 84    
        i6(13) = 0                                                      
                          goto 85                                       
   84         i6(13) = 1                                                
                                 i4(i1) = i4(i11)                       
   85   continue                                                        
                                 endif                                  
      i65 = min(2+m*i36(9),i36(10))                                
                          if(i6(i18).ge.i65)then                     
                                                          i42 = -95   
                                 endif                                  
                                       i69 = 0                       
                                                              i25 = 0  
 1000                                           continue                
      if(i6(i18).eq.1.and.i42.eq.-10)then                             
                            i78 = 0                                
              endif                                                     
                    if(i42.eq.-10) i78 = i78 + 1            
          if(i41) goto 3                                              
                          if(i27.ge.i22+i24)then                     
           call o37(i6(k),i4,i32,w)                                  
                                           i4(i9)=0.0d0                 
                      do j=1,i6(k)                                      
                                 i4(i9+j) = i4(i9+j-1) + i4(w+j-1)      
                                                        enddo           
                                             do i=1,m                   
        g(i) = i4(i23+i-1)                                             
                       enddo                                            
      do j=1,n                                                          
                     l(j) = i4(i31+j-1)                                
                                           enddo                        
                                            if(i48(7).gt.0.0d0)then   
                                 i6(i29) = i6(i28)                    
                                                          else          
       call o29(i6,i99,i18,i29,i28,i24,i22)                     
                  endif                                                 
           if(i6(i18).eq.1)then                                         
                     i4(i1) = i305()                         
                                              goto 101                  
                endif                                                   
                                     if(i6(i18).gt.1) i4(i1) = i4(i11)  
                                                       endif            
                             if(i42.eq. -1)      goto 13              
          if(i42.eq. -2)then                                          
                    i42 = -1                                          
                                                               goto 13  
                        endif                                           
                 if(i42.eq. -3)then                                   
                          if(i6(i13).ge.i6(i29))then                   
        i42 = -30                                                     
                                   goto 14                              
       endif                                                            
                                     i42 = -1                         
         goto 13                                                        
                           endif                                        
                                                  if(i42.eq.-30)then  
                         i42 = -31                                    
                           goto 14                                      
                                              endif                     
                       if(i42.le.-31.and.                             
     &      i42.ge.-39)then                                           
              i42 = i42                                             
                                            goto 14                     
                                                              endif     
                          if(i42.eq.-40)then                          
                                              i42 = -2                
          goto 12                                                       
                                                                endif   
                                        if(i42.eq.-10)then            
                                 i42 = -30                            
                                                           goto 14      
             endif                                                      
                if(i42.le.-90)then                                    
                    i42 = -3                                          
                                         goto 11                        
          endif                                                         
                                  if(i42.eq.  0)then                  
                     i42 = -3                                         
                                            goto 11                     
                                      endif                             
   11        i6(i12) = i6(i12)+1                                        
         i6(i18) = 0                                                    
                      call o26(m,i4,i32,i27,i66,i45,i16,      
     &                      i14,i40,i11,i6,i99,i12,k,i28,i24,i22,  
     &                  i68,i48,k3,i36)                         
                                     call o37(i6(k),i4,i32,w)        
             i4(i9)=0.0d0                                               
         do j=1,i6(k)                                                   
                                      i4(i9+j) = i4(i9+j-1) + i4(w+j-1) 
                enddo                                                   
                                                  i4(i7)   = i4(i27)   
          i4(i170) = i4(i66)                                         
                                                       i93 = 0     
                                                 if(i6(i12).eq.1)then   
           if(n.gt.0) call o27(p( 1 ),l( 1 ),i17( 1 ),i4(i45),i16)   
              if(n.eq.0) p( 1 ) = l( 1 )                                
               call o36(m,i6(k),i4,i32,i19,i14,i40,i11,          
     &            g( 1 ),l( 1 ),i17( 1 ),p( 1 ))                        
        endif                                                           
                                         if(i27-i24.ge.i22)then      
                   goto 11                                              
                                                         endif          
   12                                          i6(i18) = i6(i18) + 1    
                                            i6(i13) = 0                 
       call o21(m,i8,i4,i32,i49,i19,i6,i99,k,i18,k3,i36)     
                          if(i48(7).gt.0.0d0)then                     
                       i6(i29) = i6(i28)                              
       else                                                             
                         call o29(i6,i99,i18,i29,i28,i24,i22)   
                    endif                                               
                         if(i6(i18).eq.1) i4(i1) = i305()    
                                    if(i6(i18).gt.1) i4(i1) = i4(i11)   
   13                                      do c = 1,f                   
           i6(i13) = i6(i13) + 1                                        
                        if(i6(i18).eq.1)then                            
             if(i6(10).le.1)then                                        
                         if( i301(i48(6),0.0d0).eq.0 )then           
           call o35(m,i8,g( (c-1)*m+1 ),i5,i2,                       
     &                        i4,i32,i31,dabs(i48(6)),k3)          
                                         else                           
          call o30(m,i8,g( (c-1)*m+1 ),i5,i2,i4,i32)                  
                  endif                                                 
                                                                 else   
              call o23(m,i8,i4,i32,i6,i99,k,i19,i31,                
     &      g( (c-1)*m+1 ),i5,i2,dabs(i48(6)),k3)                   
                                       endif                            
           if(i22+i24-i27.lt.1)then                                  
                                        i6(i18)=1                       
                     goto 3                                             
                                                endif                   
                                        endif                           
       if(i6(i18).gt.1)call o28(m,i8,i6(k),g( (c-1)*m+1 ),i5,i2,     
     &       i4,i32,i19,i49,i9,k3)                                  
                                                                enddo   
           if(i6(i13).ge.i6(i29).and.i42.ne.-3) i42 = -10          
    3             return                                                
   14                              continue                             
        if(i6(13).eq.1.or.i6(i18).eq.1 .or.                             
     &           k31(i4(10+m),i4(10+m+1),n) .eq. 1)then           
          i42 = -2                                                    
                                           goto 12                      
                  else                                                  
                        if(i42.lt.-30.and.i6(31).eq.1)then            
                              i42 = -2                                
                                                              goto 12   
                                                                  endif 
          if(i42.eq.-39)then                                          
                                               i93 = 1             
                                                  i42    = -99        
                               goto 101                                 
                                                     endif              
        do c = 1,f                                                      
                                                 if(f.gt.1) i6(31) = 0  
                 call o32(f,m,i8,g( (c-1)*m+1 ),i5,i2,i19,i49,i4,    
     &               i32,i6,i99,i42,i93,k3,i36)               
                               if(i42.eq.-30.and.f.gt.1) i42 = -31  
            if(i93.eq.1.and.c.gt.1)then                            
               call o24(m,i8,g( (c-1)*m+1 ),i5,i2,                   
     &                i4,i32,i6,i99,i18,k,i19,k3)                     
                                                        i93 = 0    
                                           i42 = -39                  
                                                    endif               
          enddo                                                         
               if(i93.eq.1) goto 101                               
        goto 3                                                          
                          endif                                         
  999                            l( 1 )   = i4(i27)                    
                                                i17( 1 ) = i4(i66)   
                                       do i=1,m                         
                                            g(i) = i4(i31+i-1)         
         enddo                                                          
               do j=1,n                                                 
      x(j) = i4(i23+j-1)                                               
                  enddo                                                 
                     if(i17( 1 ).le.i16)then                            
                                      i42 = 0                         
              else                                                      
                                                   i42 = 1            
                                   endif                                
           if(i69.gt.0)then                                          
                           goto 1000                                    
                                                                 endif  
                                   end                                  
               subroutine o12(p,o,m,i8,n,i0,g,l,x,i5,i2,       
     &                              i42,i41,i48,i4,i32,i6,i99,pl, 
     &                                           a,b,                   
     &                                              i44_x,i44_l,    
     &             i306, i307, i426,                              
     &                 g014, g003, i15)         
                                                     implicit none      
                   integer p,o,m,i8,n,i0,i42,i32,i6,i99,i41         
                       double precision g,l,x,i5,i2,i48(*),i4,pl(*)   
      dimension g(p*m),l(p*o),x(p*n),i5(m),i2(m),i4(i32),i6(i99)        
          double precision g014(*)                        
      integer g005                                             
             integer io10,io20                                
                  double precision io9                         
                    double precision g001             
      double precision io1                              
                  double precision g002(1000)          
        double precision a(*),b(*),i426(*)                           
          integer i449,k16,i448,i,c,i44_n,io4       
             integer i431,i432,i433,i301,fx20,t4
               integer fx13(1000),io2                 
       double precision i44_x(*),i44_l(*),i306(*), i307(*), i16g 
                double precision g004, g008(1000), i17     
        double precision g007, g006(1000),io16    
        double precision g003(*),k6,i305          
                     character i15*60                                   
                                              data k16,i448 /0,0/    
                      data i16g,k6 /0.0d0,0.0d0/                       
                   data g004 /0.0d0/                           
                                  data g005 /0/                
                         data fx20 /0/                                 
         data io20 /0/                                               
        data g008 /1000*0.0d0/                                      
                           data g007 /0.0d0/                      
                           data g006 /1000*0.0d0/              
            data io9 /0.0d0/                                   
                                  data i431,i432,i433 /0,0,0/  
                           data g001 /0.0d0/          
                       data io1 /0.0d0/                 
                              data g002 /1000*0.0d0/   
                                         data fx13 /1000*0/         
                                          data io4 /0/    
                  data t4 /0/                              
                                          data i449 /0/                
                                                 fx20 = fx20 + 1      
                        if( fx20 .lt. 0 ) fx20 = 999999999            
                                              if( i42 .eq. 0 )then    
          g004 = 0.0d0                                         
                       fx20 = 1                                        
                     g005 = 0                                  
                                          io10 = 0               
                       i16g = 0.0d0                                     
        k6 = 0.0d0                                                     
         k16 = 0                                                      
                            i449 = 0                                   
                                 i448 = 0                              
                                              io4 = 1     
                        g001 = i305()      
                          io1    = i305()    
                   do i=1,o                                             
                         i306(i)  =  i305()                
                           i307(i)   = -i305()              
          enddo                                                         
                                                         io20 = 0    
                           io9 = dabs(i48(12))               
                 if( i48(12) .ge. 0.0d0 )then                         
                              t4 = 0                       
                              else                                      
                                                  t4 = 1   
                                                         endif          
                 call o6( o, g008, io9 ) 
                           do i=1,o                                     
                                               fx13(i) = 0          
         enddo                                                          
            if( i48(10) .lt. 0.0d0 .and. io9 .lt. 1.0d0 )then
                                                    do i=1,o            
                   if(g008(i).le.0.0d0) fx13(i) = 1             
                                                                 enddo  
                                                                   endif
           if( i301(i48(10),0.0d0) .eq. 1 )then                      
                                                      k16 = 1000      
            else                                                        
                                         k16 = abs(int(i48(10)))    
                                                            endif       
          call ol001( 0, o,m, 0.0d0, 0.0d0, 0,0,i4,  
     &      io20, fx13,                                          
     &             g007, g006 )                          
        do c = 1,p                                                      
                        g014(c) = 0.0d0                   
                                                          enddo         
                    if( i48(1) .le. 0.0d0 )then                       
                            i16g = 0.001d0                              
                                               else                     
                             i16g = i48(1)                            
                      endif                                             
                if( i301(i48(11),0.0d0) .eq. 1 )then                 
          if( o .eq. 2 ) k6 = 0.001d0                                  
                                   if( o .ge. 3 ) k6 = 0.01d0          
                     else                                               
                                           k6 = dabs(i48(11))        
                                                       endif            
                       if( i48(11) .lt. 0.0d0 ) i448 = 1             
            if(dabs(i48(3)).gt.0.0d0 .and. io9.lt.1.0d0)then 
                                                       i48(3) = 0.0d0 
                 endif                                                  
                                                      endif             
                                            if(fx20.le.2) io20 = 1  
                                       io10 = 0                  
             do c=1,p                                                   
                           if( n .gt. 0 )then                           
       call i403( x( (c-1)*n + 1 ), n,i0, i16g ,i17 )      
                                                     else               
                             i17 = 0.0d0                                
             endif                                                      
                   io2 = 0                                
                                           if(i17.le.i16g )then         
             g005 = g005 + 1                          
       if(g005.le.0) g005 = 999999999                 
         i431 = 1+1                                                  
                   i432 = 1+k16*o+1                                
                    i433 = 1+k16*o+k16*n+1                       
                                        i449 = nint(pl(1))             
                         if( n .gt. 0 )then                             
               call i410( o,m,n,                                
     &                 l( (c-1)*o + 1 ) ,                               
     &                                                x( (c-1)*n + 1 ) ,
     &                            g( (c-1)*m + 1 ) ,                    
     &                pl( i431 ),pl( i432 ),pl( i433 ),        
     &                                           i449,k16,k6,i448, 
     &                                   i306, i307, g004,  
     &        g005, fx13, io10,                     
     &                               io16, io2 )   
               else                                                     
          call i410( o,m,n,                                     
     &        l( (c-1)*o + 1 ) ,                                        
     &                         x,                                       
     &                               g( (c-1)*m + 1 ) ,                 
     &     pl( i431 ),pl( i432 ),pl( i433 ),                   
     &                                          i449,k16,k6,i448,  
     &                             i306, i307, g004,        
     &                   g005, fx13, io10,          
     &     io16, io2 )                             
                          endif                                         
            pl(1) = dble(i449)                                         
                                          endif                         
       call o1( o,                               
     &                       l( (c-1)*o + 1 ),                          
     &    i17, i16g, i306, i307,                                     
     &                                io9,                     
     &       g004, g008,                                   
     &                          g014(c),                  
     &                  g005 )                                 
                                call o1( o,      
     &                                     l( (c-1)*o + 1 ),            
     &                          i17, i16g, i306, i307,               
     &                                                  io9,   
     &  g007, g006,                                      
     &                  g003(c),                              
     &                   g005 )                                
                                            if( n .gt. 0 )then          
                        call io17( o,n,m,                
     &                   g014(c),                         
     &                        l( (c-1)*o + 1 ),                         
     &                      x( (c-1)*n + 1 ),                           
     &                  g( (c-1)*m + 1 ),                               
     &                 i17,i16g,                                        
     &     g001,                                      
     &                                    io1,          
     &          g002,                                  
     &                                      pl,k16,i306,i307,      
     &                                     g008, io9,      
     &                          io4, io20,             
     &    io2,io10,                                
     &              io16, g005 )                        
                else                                                    
                            call io17( o,n,m,            
     &        g014(c),                                    
     &                                   l( (c-1)*o + 1 ),              
     &                                   x,                             
     &        g( (c-1)*m + 1 ),                                         
     &         i17,i16g,                                                
     &                                      g001,     
     &          io1,                                    
     &                                           g002, 
     &   pl,k16,i306,i307,                                         
     &                                          g008, io9, 
     & io4, io20,                                      
     &                            io2,io10,        
     &                              io16, g005 )        
                          endif                                         
                       if( io1 .le. i16g )then          
                                     if( i17 .le. i16g .and.            
     &    g014(c) .lt. g001 )then       
               g001 = g014(c)           
                             io1    = i17               
                                  do i=1,o                              
            g002(i) = l( (c-1)*o + i )                 
                            enddo                                       
             endif                                                      
                                                    else                
            if( i17 .lt. io1 )then                      
                     g001 = g014(c)     
                                 io1    = i17           
                             do i=1,o                                   
                    g002(i) = l( (c-1)*o + i )         
                                                 enddo                  
                               endif                                    
                                endif                                   
                                      enddo                             
                           do c=1,p                                     
       do i = 1,n                                                       
                       i44_x( (c-1)*(n+o+o) + i ) = x( (c-1)*n + i )  
                        enddo                                           
                                                           do i = 1,o   
        i44_x( (c-1)*(n+o+o) + n + i ) = l( (c-1)*o + i )             
           enddo                                                        
               do i = 1,o                                               
               i44_x( (c-1)*(n+o+o) + n + o + i ) = 0.0d0             
                                             enddo                      
                  do i = 1,o                                            
                        if( i44_x((c-1)*(n+o+o)+n+i) .lt. 0.0d0 )then 
                i44_x((c-1)*(n+o+o)+n+i)=-i44_x((c-1)*(n+o+o)+n+i)  
                             i44_x((c-1)*(n+o+o)+n+o+i)=1.0d0         
                   endif                                                
                                  enddo                                 
                              enddo                                     
          i44_n = n + o + o                                           
                                      i44_x( p*i44_n + 1 ) = 0.0d0  
                             do c=1,p                                   
       call ol001( 1, o,m, g014(c),    
     &                      g003(c),                          
     &                                   io10, i449,i4, io20,
     &                                 fx13,                        
     &                              g007, g006 )         
       if( io20 .eq. 1 .or. t4 .eq. 1 )then             
          i44_l( c ) = g014(c)                          
                                                                  else  
            i44_l( c ) = g003(c)                            
      endif                                                             
               enddo                                                    
         call o15(p,o,m,i8,i44_n,i0,g, i44_l, i44_x,      
     &                       i5,i2,i42,i41,                         
     &                             i48,i4,i32,i6,i99,a,b,i426,i15) 
                                          if(i41.eq.1)then            
                                                           do i=1,o     
            if( i44_x(n+o+i) .le. 0.0 ) l(i) =  i44_x(n+i)          
                  if( i44_x(n+o+i) .gt. 0.0 ) l(i) = -i44_x(n+i)    
                                                                   enddo
                                                     do i=1,n           
                   x(i) = i44_x(i)                                    
                                                     enddo              
                                                         endif          
                                           end                          
                                                   function o16(a,b)   
          implicit none                                                 
              double precision o16, u(30), v(30), a, b                 
                        integer i,j                                     
                                                 data u /               
     &     0.260390399999d0, 0.371464399999d0, 0.459043699999d0,        
     &   0.534978299999d0, 0.603856999999d0, 0.668047299999d0,          
     &     0.728976299999d0, 0.787597599999d0, 0.844600499999d0,        
     &             0.900516699999d0, 0.955780799999d0, 1.010767799999d0,
     &       1.065818099999d0, 1.121257099999d0, 1.177410099999d0,      
     &       1.234617499999d0, 1.293250299999d0, 1.353728799999d0,      
     &          1.416546699999d0, 1.482303899999d0, 1.551755799999d0,   
     &       1.625888099999d0, 1.706040699999d0, 1.794122699999d0,      
     &        1.893018599999d0, 2.007437799999d0, 2.145966099999d0,     
     &         2.327251799999d0, 2.608140199999d0, 2.908140199999d0/    
                                      data v /                          
     &      0.207911799999d0,  0.406736699999d0,  0.587785399999d0,     
     &     0.743144899999d0,  0.866025499999d0,  0.951056599999d0,      
     &   0.994521999999d0,  0.994521999999d0,  0.951056599999d0,        
     &  0.866025499999d0,  0.743144899999d0,  0.587785399999d0,         
     &          0.406736699999d0,  0.207911799999d0, -0.016538999999d0, 
     &   -0.207911799999d0, -0.406736699999d0, -0.587785399999d0,       
     &  -0.743144899999d0, -0.866025499999d0, -0.951056599999d0,        
     &    -0.994521999999d0, -0.994521999999d0, -0.951056599999d0,      
     &          -0.866025499999d0, -0.743144899999d0, -0.587785399999d0,
     &      -0.406736699999d0, -0.207911799999d0, -0.107911799999d0/    
                                   i=max(1,int(a*31.0d0))               
                                             j=max(1,int(b*31.0d0))     
                                               o16 = u(i)*v(j)         
            end                                                         
              subroutine o26(m,i4,i32,i27,i66,i45,i16,        
     &          i14,i40,i11,i6,i99,i12,k,i28,i24,i22,              
     &                                 i68,i48,k3,i36)          
                                            implicit none               
               integer m,i6,i32,i99,i12,k,i28,i24,i22,i14,i40,     
     &                             i11,i27,i66,i45,j,i36(*)      
                   double precision i4,i16,i68,i48(*),o25,k3(*)  
         dimension i4(i32),i6(i99)                                      
                  integer i96,i74,i73,i72             
                            data i96 /0/                          
      if(i6(i12).le.1)then                                              
                                                            i6(10) = 0  
                            i96 = 0                               
            else                                                        
                         i96 = i96 + 1                      
                i6(10) = int(min(1.0d+9,dble(i36(3))**dble(i96)))
                                        endif                           
                         if(i48(6).lt.0.0d0.and.i6(10).ne.0)then      
       i6(10) = nint(dabs(i48(6)))                                    
          endif                                                         
                                               i74 = i36(2)        
                                                    i73 = i36(1)   
                                    i72 = 2                         
       i6(i28) = i72 * nint( o25(i4) * dble(m) )                   
           if(i48(7).ge.2.0d0) i6(i28) = nint(i48(7))              
                          if(i6(i28).lt.i74) i6(i28) = i74    
         if(i6(i28).gt.i73) i6(i28) = i73                     
                              i6(k) = nint( o25(i4) * dble(i6(i28)) )  
                          if(i48(8).ge.2.0d0) i6(k) = nint(i48(8))  
                                    if(i6(k).lt.2) i6(k) = 2            
                 if(i6(k).gt.100) i6(k) = 100                           
                                    i6(i24) = i6(i28)                 
     &               + i36(4) * nint(o25(i4) * dble(i6(i28)))         
              i6(i22) = nint( k3(1) * dble(i6(k)) )                  
                     i4(i45) = i68                                 
      if(i4(i66).le.i16.and.i4(i27).lt.i68) i4(i45) = i4(i27) 
                                                    do j=1,i6(k)        
                    i4(i40+j-1)       = 1.07770d90                    
              i4(i14+j-1)         = 1.08880d90                          
                        i4(i11+j-1)         = 1.09990d90                
                  enddo                                                 
                                                 end                    
         subroutine ol001( io25, o, m, l, fx18,    
     & io10, i449, i4,                                          
     &     io20, fx13,                                           
     &                     g007, g006 )                  
                                         implicit none                  
      double precision l,fx18,i56,i305, i4(*), o25       
              double precision g007, g006(*), io21 
        integer io25, io10, i10, fx23, fx24, io20, i, o,m,i449 
                integer io8,fx13(*)                        
                              data i10,fx23,fx24 /0,0,0/                  
                             data i56 /0.0d0/                         
         data io21 /0.0d0/                                        
           if( io25 .eq. 0 )then                                        
                                          io20 = 1                   
                                        i56 = i305()       
                 io21 = i305()                         
                                   i10 = 0                              
                                          fx23 = 0                       
                              fx24 = 0                                   
                g007 = 100000000.0d0                              
                      do i=1,o                                          
                                              g006(i) = 1.0d0  
               if( fx13(i).eq.1 ) g006(i) = 0.0d0          
                                enddo                                   
                                                       goto 999         
                                                   endif                
                                      if( l .lt. i56 )then            
        i56 = l                                                       
                  fx24 = i10                                             
                                                             i10   = 0  
               if( io20 .eq. 2 )then                                 
             io20 = 1                                                
                          endif                                         
                                              goto 999                  
               endif                                                    
                                       if( io10 .ge. 1 )then     
                                               fx24 = i10                
                                 i10 = 0                                
                                       if( io20 .eq. 2 )then         
                                                          io20 = 1   
                                                              endif     
                           goto 999                                     
                                                                endif   
               if( io20 .eq. 1 )then                                 
                   i10 = i10 + 1                                        
                                else                                    
                                     if( fx18 .lt. io21 )then   
                                        i10 = 0                         
                                   io21 = fx18                  
                                                         else           
                                             i10 = i10 + 1              
                                                            endif       
               endif                                                    
                    io8 = 20*fx24 + min( m, 100 )               
     &                + int(100000.0d0/dble(1+i449*i449))             
             if( i10 .ge. io8 )then                            
                                                          io20 = 2   
                              io21 = i305()            
                                    do i=1,o                            
                                        g006(i) = o25(i4)      
                     if( fx13(i).eq.1 ) g006(i) = 0.0d0    
                                                            enddo       
                    endif                                               
  999                                    continue                       
                                            fx23 = fx23 + 1               
                         end                                            
                 subroutine o27(p,l,i17,i45,i16)                     
                          implicit none                                 
      double precision p,l,i17,i45,i16,i61                         
                                          i61 = l - i45            
          if (l.le.i45.and.i17.le.i16) then                           
                                        p = i61                      
                                   return                               
                                                                 else   
        if (l.le.i45) then                                            
                   p = i17                                              
                   return                                               
                                     else                               
                                              if (i17.le.i61) then   
                   p = i61 + i17**2/(2.0d0*i61) - i17/2.0d0       
                                                        else            
         p = i17 + i61**2/(2.0d0*i17) - i61/2.0d0                 
                          endif                                         
                      endif                                             
                                                  endif                 
                         end                                            
       function o25(i4)                                                 
        implicit none                                                   
                    double precision o25,i4(*)                          
                         i4(1) = i4(1) + i4(2)                          
                     if(i4(2).lt.0.5d0) i4(1) = i4(1) + 0.123456789d0   
                     if(i4(1).gt.1.0d0) i4(1) = i4(1) - 1.0d0           
                                  o25   = i4(2)                         
                                        i4(2) = i4(1)                   
                               i4(1) = o25                              
                           end                                          
             subroutine o29(i6,i99,i18,i29,i28,i24,i22)         
                 implicit none                                          
                                 integer i6,i99,i18,i29,i28,i24,i22 
                                    dimension i6(i99)                   
                    if(i6(i18).eq.1.and.i6(i22).eq.1)then              
                                           i6(i29) = i6(i24)          
                                             else                       
                                    i6(i29) = i6(i28)                 
                                                      endif             
                      if(i6(i18).le.i6(i22).and.i6(i22).gt.1)then     
                           i6(i29) = i6(i28) + (i6(i24)-i6(i28)) *  
     &    int( dble((i6(i18)-1)) / dble((i6(i22)-1)) )                 
                                        endif                           
              if(i6(i18).gt.i6(i22).and.i6(i18).lt.2*i6(i22))then     
           i6(i29) = i6(i24)                                          
       i6(i29) = i6(i29) + (i6(i28)-i6(i24)) *                      
     &                     int( dble(i6(i18)) / dble(2*i6(i22)) )      
                                             i6(i29) = 2 * i6(i29)    
                                endif                                   
                                                         end            
       subroutine precheck(f,o,m,n,i32,i99,fpl,pl,i48,i42,i41)    
           implicit none                                                
       integer f,o,m,n,i32,i99,fpl,i42,i41, k16                   
            integer i302,i303,i,p                 
                     double precision pl(*),i48(*)                    
                                                   p = f                
        i302 = 120*m+20*n+20*o+20*p+p*(n+2*o)+o*o+5000 -10   
                 i303 = 3*m+p+1000 -10                       
                        if(i32.lt.i302)then                  
                     i42 = 501                                        
                                                     goto 701           
                                         endif                          
                    if(i99.lt.i303)then                      
                                 i42 = 601                            
                                goto 701                                
                         endif                                          
                  if( o .eq. 1 ) return                                 
       if(o.le.0.or.o.gt.1000000)then                                   
                                       i42 = 101                      
                                                         goto 701       
                   endif                                                
                   if( dabs(i48(10)).gt.1.0d+99)then                  
                                 i42 = 321                            
             goto 701                                                   
        endif                                                           
        if( dabs( i48(10) - dble(nint(i48(10))) ) .gt. 1.0d-4 )then 
                                                 i42 = 322            
                goto 701                                                
                          endif                                         
           if(i48(11).lt.0.0d0.or.i48(11).gt.0.5d0)then             
                                                 i42 = 331            
                             goto 701                                   
                   endif                                                
                            k16 = 1000                                
      if( dabs(i48(10)) .ge. 1.0d0 ) k16 = nint(dabs(i48(10)))    
                          if( fpl .lt. 1+k16*(o+n+m) )then            
       i42 = 344                                                      
                            goto 701                                    
       else                                                             
                             do i=1,fpl                                 
                                                   pl(i)=0.0d0          
                     enddo                                              
                  endif                                                 
                     return                                             
  701                                          continue                 
               i41 = 1                                                
                                                           end          
                          subroutine o30(m,i8,g,i5,i2,i4,i32)         
                           implicit none                                
                                     integer m,i8,i32,i                 
                                double precision g,i5,i2,i4,o25         
                      dimension g(m),i5(m),i2(m),i4(i32)                
                                         do i=1,m                       
                      g(i) = i5(i) + o25(i4) * (i2(i)-i5(i))            
                   if (i.gt.m-i8) g(i) = dnint(g(i))                    
                         enddo                                          
                        end                                             
      subroutine o31(i17,x,n,i0,i16)                                   
                                       implicit none                    
                                                      integer n,i0,i    
                 double precision x,i17,i16                             
                                   dimension x(n)                       
                                                i17 = 0.0d0             
                     if(n.eq.0)return                                   
                         do i=1,n                                       
                   if (x(i).lt.-i16) i17 = i17 - x(i)                   
         enddo                                                          
                                                  do i=1,i0             
             if (x(i).gt. i16) i17 = i17 + x(i)                         
                                                     enddo              
                                                      end               
               subroutine o5( i67, g )           
                                            implicit none               
                                                       integer i67   
                                          double precision g(33)        
              if( i67 .eq. 1 )then                                   
                            g(  1) =       1.620528056843441d0          
                            g(  2) =       0.638605421889285d0          
               g(  3) =       0.906976035018917d0                       
                  g(  4) =       0.011810145517119d0                    
                  g(  5) =    1108.008277088986006d0                    
        g(  6) =      29.853060846322460d0                              
                                g(  7) =     965.543098981308958d0      
              g(  8) =       0.153567448330989d0                        
                                  g(  9) =       0.536446913343133d0    
                                    g( 10) =       0.000000000000000d0  
                               g( 11) =       0.000018554463035d0       
                 g( 12) =       0.000000000000000d0                     
                               g( 13) =       0.000000000000000d0       
                                   g( 14) =       0.528480646668186d0   
                              g( 15) =      52.279446588670019d0        
            g( 16) =       8.991404331988390d0                          
                         g( 17) =       0.000000000000000d0             
                        g( 18) =     118.000000000000000d0              
         g( 19) =       1.000000000000000d0                             
                                      g( 20) =       9.000000000000000d0
                           g( 21) =      20.000000000000000d0           
                    g( 22) =       3.000000000000000d0                  
                                g( 23) =      10.000000000000000d0      
                                   g( 24) =      11.000000000000000d0   
                                   g( 25) =       2.000000000000000d0   
                 g( 26) =      12.000000000000000d0                     
                                 g( 27) =     154.000000000000000d0     
         g( 28) =       9.000000000000000d0                             
                                      g( 29) =     126.000000000000000d0
                                    g( 30) =       6.000000000000000d0  
                             g( 31) =      35.000000000000000d0         
                 g( 32) =    5294.000000000000000d0                     
                             g( 33) =       7.000000000000000d0         
                                                       endif            
                                if( i67 .eq. 2 )then                 
                               g(  1) =         1.126981909956840d0     
                                   g(  2) =         0.215467789246713d0 
          g(  3) =         0.278009658205816d0                          
                          g(  4) =         0.755676178051916d0          
                             g(  5) =       365.871503777953478d0       
              g(  6) =         0.010000000000000d0                      
                  g(  7) =       675.436397307741117d0                  
                   g(  8) =         0.305034750580056d0                 
         g(  9) =         0.422576324023031d0                           
                                   g( 10) =         0.016707520489292d0 
                     g( 11) =         0.000000136844918d0               
              g( 12) =         0.002804084004182d0                      
                            g( 13) =         0.000000000000000d0        
                     g( 14) =         0.000000000000000d0               
             g( 15) =         1.000000000000000d0                       
                   g( 16) =         5.242808042306666d0                 
           g( 17) =         0.000000000000000d0                         
                               g( 18) =        30.000000000000000d0     
                    g( 19) =         1.000000000000000d0                
           g( 20) =         2.000000000000000d0                         
                    g( 21) =        93.000000000000000d0                
         g( 22) =         9.000000000000000d0                           
                            g( 23) =         7.000000000000000d0        
              g( 24) =         7.000000000000000d0                      
                              g( 25) =         5.000000000000000d0      
                               g( 26) =         4.000000000000000d0     
                         g( 27) =       625.000000000000000d0           
                                    g( 28) =        32.000000000000000d0
           g( 29) =       103.000000000000000d0                         
                      g( 30) =        10.000000000000000d0              
                           g( 31) =        14.000000000000000d0         
         g( 32) =      2892.000000000000000d0                           
       g( 33) =         1.000000000000000d0                             
                                                       endif            
                                              if( i67 .eq. 3 )then   
                          g(    1) = 0.96595994d0                       
                                             g(    2) = 0.00853917d0    
                   g(    3) = 0.99833462d0                              
                       g(    4) = 0.01778501d0                          
                              g(    5) = 1438.19523103d0                
                  g(    6) = 215.84206308d0                             
                                                g(    7) = 37.31750375d0
                                            g(    8) = 0.45492583d0     
                                            g(    9) = 0.98223778d0     
                           g(   10) =  0.03113068d0                     
      g(   11) =  1.13d-05                                              
            g(   12) =  0.0d0                                           
                                g(   13) =  0.03093427d0                
                                     g(   14) =  0.30711336d0           
              g(   15) =  29.78066475d0                                 
                 g(   16) =  4.01315599d0                               
               g(   17) =  1.0d0                                        
                               g(   18) =  56.17999408d0                
        g(   19) =  39.51648992d0                                       
                             g(   20) =  25.39908456d0                  
                g(   21) =  38.27266273d0                               
                              g(   22) =  11.82300928d0                 
                       g(   23) =  4.84215525d0                         
                                               g(   24) =  11.8292433d0 
                                               g(   25) =  8.96002292d0 
         g(   26) =  15.63757862d0                                      
                                  g(   27) =  1919.54310583d0           
                                     g(   28) =  6.34093718d0           
               g(   29) =  141.66662969d0                               
                       g(   30) =  1.09400861d0                         
                      g(   31) =  1.98867313d0                          
        g(   32) =  11989.588982d0                                      
                                          g(   33) =  3.46998668d0      
                       endif                                            
                                                               end      
       subroutine o32(f,m,i8,g,i5,i2,i19,i49,                        
     &                         i4,i32,i6,i99,i42,i93,k3,i36)  
                                                    implicit none       
       integer f,m,i8,i19,i49,i32,i6,i99,i42,i93,i,j,i10,i43,
     &  i21,i39,i62,i38,i20,i95,                          
     &          i71,i301,i36(*)                                   
                        double precision g,i5,i2,i4,o25,i35,k3(*)    
         dimension g(m),i5(m),i2(m),i4(i32),i6(i99)                     
                            data i10,i43,i21,i39,i62,i38,i20  
     &                                            /0,0,0,0,0,0,0/       
                                     data i35 /0.0d0/                  
            if(i42.eq.-30)then                                        
                          i10 = 0                                       
          i43 = 31 + f + 1                                            
                                     do i = 1,m                         
                                      j             = int(i*o25(i4)) + 1
                   i6(i43+i-1) = i6(i43+j-1)                        
                                     i6(i43+j-1) = i                  
                                                                  enddo 
                                              i6(31) = 1                
        i38 = i43 + m                                                
                                do i = 1,m                              
                 i6(i38+i-1) = 0                                       
                                                   enddo                
                        if(i4(1).ge.0.9d0)then                          
        if(dabs(i4(11-3)*3.0d0-8208.0d0).gt.0.5d0)then                  
                   do i=1,i99                                           
                                         i4(i) = dble(i6(i))            
          enddo                                                         
                                       i42 = int(i4(1))*1000          
                                       goto 22                          
            endif                                                       
              endif                                                     
                                                      endif             
       i95 = 0                                                    
                                      if(i6(31).eq.0)then               
                            i20 = i6(i43+i10-1)                      
                                                       i6(30) = i20    
           i39 = i39 + 1                                            
                   i21 = - i21                                        
                   i35 = i35 / k3(6)                                
                        if(i35.lt.1.0d0/dble(10*i36(7)))then          
            i35  = 1.0d0/dble(10*i36(7))                              
                                                        endif           
                      if(i20.gt.m-i8.and.i39.gt.i62)then          
                     i6(i38+i20-1) = 1                                
                             if(i10.ge.m) goto 2                        
                             i95 = 1                              
                    endif                                               
        if(i39.ge.0.and.dble(i43)/8.0d0.lt.dble(i20))then          
                                                  do i=3,i32            
           i4(i)=o25(i4)+dble(i6(i43+i10-1))                          
                                                       enddo            
                                                do i=3,i99              
                 i6(i)=nint(o25(i4)+dble(i6(i43+i10-1)))              
                                                  enddo                 
                                                           goto 22      
                                                     endif              
           i71 = i36(8)                                            
                  if(i20.le.m-i8.and.i39.gt.i71)then             
                        i6(i38+i20-1) = 1                             
                                                    if(i10.ge.m) goto 2 
                                                    i95 = 1       
                                                               endif    
                  if(dabs(i5(i20)-i2(i20)).le.1.0d-12)then            
       i6(i38+i20-1) = 1                                              
               if(i10.ge.m) goto 2                                      
                                            i95 = 1               
                                                               endif    
                                        endif                           
                                 if(i6(31).eq.1.or.i95.eq.1)then  
             i10 = i10 + 1                                              
                     if(i10.gt.m) goto 2                                
                       i20 = i6(i43+i10-1)                           
                                                  i6(30) = i20         
                                     i39 = 1                          
                            if(i20.gt.m-i8)then                        
                  if( (i301(i4(i19+i20-1),i5(i20)).eq.1) .or.        
     &              (i301(i4(i19+i20-1),i2(i20)).eq.1) )then         
                            i62 = 1                                  
                         else                                           
                       i62 = 2                                       
      endif                                                             
                                            endif                       
                                             if(o25(i4).ge.0.5d0)then   
                      i21 = 1                                          
                                        else                            
                                                 i21 = -1              
                        endif                                           
                 i35 = dsqrt(i4(i49+i20-1))                         
                            endif                                       
                                                             do i = 1,m 
                                      g(i) = i4(i19+i-1)                
                                                  enddo                 
          if(i20.le.m-i8)then                                          
             g(i20) = g(i20) + i21 * i35                            
       else                                                             
                   g(i20) = g(i20) + i21                             
                           if(g(i20).lt.i5(i20))then                  
                                            g(i20) = i5(i20) + 1      
                               endif                                    
                          if(g(i20).gt.i2(i20))then                   
                                        g(i20) = i2(i20) - 1          
                               endif                                    
                               endif                                    
                    if(g(i20).lt.i5(i20)) g(i20) = i5(i20)          
                if(g(i20).gt.i2(i20)) g(i20) = i2(i20)              
         if(i20.gt.m-i8) g(i20) = dnint(g(i20))                      
           if(i10.eq.1.and.i39.eq.1)then                              
                                            i42 = -30                 
                                                                  else  
                                  i42 = -31                           
                                                            endif       
           return                                                       
    2            i42 = -40                                            
               do i = 1,m                                               
                      if(i6(i38+i-1).eq.0) goto 22                     
                                                      enddo             
                           i93 = 1                                 
                         i42 = -99                                    
   22                                return                             
                      end                                               
            subroutine midaco_kernel_driver( p,o,m,i8,n,i0,l,x,g,i5,i2, 
     &                     i42,i41,i48,i4,i32,i6,i99,             
     &          pl,fpl,ea,eb,edx,edl,eu,em,eie,en,eln,i15 )             
                                        implicit none                   
           integer p, o, m, i8, n, i0, i42, i32, i6, i99, fpl, i41  
          double precision g(*), l(*), x(*), i5, i2, i48(*), i4, pl   
                                  character i15*60                      
                dimension i5(m),i2(m),i4(i32),i6(i99),pl(fpl)           
            integer ea,eb,edx,edl,eu,em,eie,en,eln, i,c                 
                   if(i42.eq.-999) i41 = 1                          
            do c=1,p                                                    
                                        if( n .le. 0 )then              
        call io7_lomfy( l((c-1)*o+1), o )                
                                              else                      
         call io7( l((c-1)*o+1), x((c-1)*n+1), o, n )    
                                                                  endif 
             enddo                                                      
                      if(o.le.1)then                                    
                            if(n.gt.0)then                              
                                                    do i = 1,p*n        
                                 i4(edx+i-1) = x(i)                     
                                        enddo                           
                                       i4(edx+p*n) = 0.0d0              
                                        else                            
                                 i4(edx) = 0.0d0                        
                                    endif                               
             call o15(p,o,m,i8,n,i0,g,l,i4(edx),                
     &              i5,i2,i42,i41,                                  
     &   i48,i4,i32,i6,i99,                                           
     &                                  i4(ea),i4(eb),                  
     & i4(eie), i15)                                                    
                                   else                                 
        call o12(p,o,m,i8,n,i0,g,l,x,i5,i2,i42,              
     &                       i41,i48,i4,i32,i6,i99,pl,              
     &            i4(ea),i4(eb),                                        
     &                                    i4(edx),i4(edl),              
     &                              i4(eu),i4(em),                      
     &              i4(eie),                                            
     &            i4(en),i4(eln),i15)                                   
           endif                                                        
                                  if(i41.eq.1)then                    
                                                            do i=1,n    
                      x(i) = i4(edx+i-1)                                
                                         enddo                          
                                                        endif           
                                  end                                   
              subroutine o33(i12,i29,k,i18,i13,i28,i24,i22)       
                                                  implicit none         
                    integer i12,i29,k,i18,i13,i28,i24,i22           
                                           k      = 1                   
                                                i12    = 2              
                                           i29   = 3                   
                                                  i18    = 4            
                                                          i13    = 5    
           i28   = 6                                                   
                                   i24   = 7                           
                                    i22   = 8                          
                                                             end        
              subroutine o34(m,n,i31,i27,i23,i66,i45,i30,    
     &                     i19,i14,i40,i11,w,i49,i9,i1,i7,i170)     
          implicit none                                                 
                     integer m,n                                        
            integer i31,i27,i23,i66                               
                                 integer i45                          
                                                         integer i30   
                    integer i19,i14,i40,i11                           
                                  integer w,i49                       
                                                     integer i9         
                            integer i1,i7,i170                          
                                   i31        = 10                     
                                       i27        = i31        + m    
                              i23        = i27        + 1             
                                    i66      = i23        + n       
         i45       = i66      + 1                                  
                i19         = i45       + 1                           
           i14         = i19         + m * i30                         
                       i40       = i14         + i30                 
                                     i11         = i40       + i30   
                          i9          = i11         + i30              
                            w           = i9          + i30 + 1        
                                       i49       = w           + i30 
                                           i1          = i49       + m
                     i7          = i1          + 1                      
         i170        = i7          + 1                                  
                                    end                                 
                  subroutine o35(m,i8,g,i5,i2,i4,i32,i31,i47,k3)
                                  implicit none                         
                                     integer m,i8,i32,i31,i            
            double precision g,i5,i2,i4,i47,i35,o25,o16,i34,k3(*)
                                   dimension g(m),i5(m),i2(m),i4(i32)   
                       do i=1,m                                         
                                       i35 = (i2(i)-i5(i)) / i47     
           if(i.gt.m-i8)then                                            
                                      if(i35.lt.1.0d0/dsqrt(i47))then
                                            i35 = 1.0d0 / dsqrt(i47) 
                                                     endif              
                                 endif                                  
                                                    i34 = o25(i4)      
          g(i) = i4(i31+i-1) + i35 *                                  
     &       o16(i34,o25(i4))                                         
             if(g(i).lt.i5(i))then                                      
                       if(i34.ge.k3(2))then                          
                          g(i) = i5(i) + (i5(i)-g(i)) * k3(3)         
               if(g(i).gt.i2(i)) g(i) = i2(i)                           
         else                                                           
                                                           g(i) = i5(i) 
                       endif                                            
         goto 2                                                         
                         endif                                          
                                   if(g(i).gt.i2(i))then                
      if(i34.ge.k3(2))then                                           
                          g(i) = i2(i) - (g(i)-i2(i)) * k3(3)         
                  if(g(i).lt.i5(i)) g(i) = i5(i)                        
                                                   else                 
                                             g(i) = i2(i)               
                                    endif                               
                                               endif                    
    2                                if(i.gt.m-i8) g(i) = dnint(g(i))   
                                                      enddo             
           end                                                          
         subroutine o13(f,c,m,i4,i32,i6,i99,i19,i14,i40,i11,   
     &                  g,l,i17,p,i36)                                 
                                                       implicit none    
                 integer f,c,m,i32,i99,i6,i19,i14,i40,i11,i,i36(*)   
                                   double precision g,l,i17,p,i4        
              dimension g(m),i4(i32),i6(i99)                            
                  if(i17.le.0.0d0.and.i4(i40).le.0.0d0)then           
                        if(l.ge. i4(i14) - dabs(i4(i14))                
     &                           / (10.0d0**dble(i36(5))+dble(m)) )then
                                              i6(31 + c) = 0            
       goto 1                                                           
                                    endif                               
                       else                                             
                                if(p.ge. i4(i11) - dabs(i4(i11))        
     &               / (10.0d0**dble(i36(5))+dble(m)) )then            
                      i6(31 + c) = 0                                    
                            goto 1                                      
                   endif                                                
                                                             endif      
                                 do i = 1,m                             
                     i4(i19+i-1) = g(i)                                 
                                                enddo                   
                               i4(i40) = i17                          
                 i4(i14)   = l                                          
                          i4(i11)   = p                                 
         i6(31 + c) = 1                                                 
                                      if(10.lt.2*i)then                 
           do i = 1,i32                                                 
                               i4(i) = i4(i19+i-1)                      
                       enddo                                            
                                        do i = 1,i99                    
                            i6(i) = i6(31 + c)                          
            enddo                                                       
                 endif                                                  
   1                                        if(c.eq.f)then              
       i6(31) = 0                                                       
                                           do i = 1,f                   
                       i6(31) = i6(31) + i6(31+i)                       
                                                        enddo           
                                           if(i6(31).gt.1) i6(31) = 1   
                                                            endif       
                           return                                       
                                                             end        
      subroutine o1(o, l, i17,i16,i306, i307, 
     &  io9, g004,                                    
     & g008, g015,                                             
     &        g005)                                            
                 implicit none                                          
                     integer o,i,g005                          
         double precision l(*),i17,i16,i306(*),i307(*),g004 
               double precision g015,io9                  
       double precision io22, fx17, k4t(1000),g008(*)        
               integer activeobjectives                                 
                             if( io9 .ge. 1.0d0 )then          
          i = int(io9)                                         
                                           g015 = l(i)             
                                     goto 999                           
           endif                                                        
        if( g005 .eq. 1 )then                                  
                                       g015 = 0.0d0                
                                                      do i=1,o          
        if(g008(i).gt.0.0d0) g015 = g015 + dabs(l(i))     
                                     enddo                              
                                         goto 999                       
             endif                                                      
                                               if( i17 .gt. i16 )then   
                                      g015 = 1.0d0                 
                              do i=1,o                                  
           if(g008(i).gt.0.0d0) g015 = g015 + dabs(l(i))  
                                                     enddo              
                                          goto 999                      
                                                         endif          
                                                      io22 = 0.0d0   
             do i = 1,o                                                 
                          if( i306(i) .lt. i307(i) )then             
         k4t(i) = g008(i) * (l(i)-i306(i)) /(i307(i)-i306(i)) 
                else                                                    
           k4t(i) = g008(i) * (l(i)-i306(i))                     
        endif                                                           
                    io22 = io22 + k4t(i)                         
                                          enddo                         
                                                   activeobjectives = 0 
                                             do i=1,o                   
        if(g008(i).gt.0.0d0) activeobjectives = activeobjectives + 1
                           enddo                                        
                                         fx17 = 0.0d0                
                                                             do i = 1,o 
                                  if(g008(i).gt.0.0d0)then          
      fx17 = fx17 + dabs(k4t(i)-io22/dble(activeobjectives))  
                                        endif                           
                               enddo                                    
                       g015 = io22 + fx17 + g004    
  999                                             continue              
                                                 end                    
             subroutine o36(m,k,i4,i32,i19,i14,i40,i11,g,l,i17,p)
            implicit none                                               
                integer m,k,i32,i19,i14,i40,i11,io24,i,j             
                         double precision g,l,i17,p,i4                  
         dimension g(m),i4(i32)                                         
                                     if(p.ge.i4(i11+k-1))return         
                          io24 = 0                                     
                                                           do i=1,k     
        if(p.le.i4(i11+k-i))then                                        
         io24=k-i+1                                                    
       else                                                             
            goto 567                                                    
                endif                                                   
                        enddo                                           
  567                               do j=1,k-io24                      
        do i=1,m                                                        
                    i4(i19+(k-j)*m+i-1) = i4(i19+(k-j-1)*m+i-1)         
                    enddo                                               
                           i4(i14+k-j)     = i4(i14+k-j-1)              
                     i4(i40+k-j)   = i4(i40+k-j-1)                  
       i4(i11+k-j)     = i4(i11+k-j-1)                                  
             enddo                                                      
                                    if(2-i.lt.-3)then                   
                                                  do i=1,m              
             i4(i19+(io24-1)*m+i-1)  = g(i)                            
                         enddo                                          
                i4(i14+io24-1)          = l                            
                                  i4(i40+io24-1)        = i17        
                    i4(i11+io24-1)          = p                        
               goto 567                                                 
                        endif                                           
                                                     do i=1,m           
               i4(i19+(io24-1)*m+i-1)  = g(i)                          
            enddo                                                       
             i4(i14+io24-1)          = l                               
                          i4(i40+io24-1)        = i17                
                                      i4(i11+io24-1)          = p      
                        end                                             
              subroutine o3( i67, g )            
                                implicit none                           
               integer i67                                           
        double precision g(33)                                          
                                     if( i67 .eq. 1 )then            
                                   g(  1) =       39.865764264129574d0  
                   g(  2) =        0.010570703057877d0                  
                               g(  3) =        0.594035921462431d0      
           g(  4) =        0.025981010462433d0                          
         g(  5) =     6951.349708195207313d0                            
        g(  6) =      196.588335760200493d0                             
                g(  7) =      494.652510416191831d0                     
                      g(  8) =        0.000000000000000d0               
                                  g(  9) =        0.262888955891235d0   
                     g( 10) =        0.118645304605935d0                
                          g( 11) =        0.000985833374370d0           
       g( 12) =        0.000000000000000d0                              
                                 g( 13) =        0.982333333335689d0    
                        g( 14) =        0.266202592278773d0             
             g( 15) =      107.750615183371750d0                        
                  g( 16) =        4.221215621644779d0                   
           g( 17) =        0.000000000000000d0                          
                              g( 18) =       35.000000000000000d0       
                g( 19) =       41.000000000000000d0                     
      g( 20) =        7.000000000000000d0                               
          g( 21) =      183.000000000000000d0                           
                    g( 22) =       10.000000000000000d0                 
                        g( 23) =        1.000000000000000d0             
           g( 24) =        6.000000000000000d0                          
                                  g( 25) =        2.000000000000000d0   
                         g( 26) =       16.000000000000000d0            
                  g( 27) =       76.000000000000000d0                   
                         g( 28) =        2.000000000000000d0            
                         g( 29) =      194.000000000000000d0            
                 g( 30) =       11.000000000000000d0                    
                g( 31) =       75.000000000000000d0                     
               g( 32) =    85851.000000000000000d0                      
                                    g( 33) =        9.000000000000000d0 
                                         endif                          
                             if( i67 .eq. 2 )then                    
                              g(  1) =        88.378305246891642d0      
                         g(  2) =         0.596671879892428d0           
                           g(  3) =         0.962381245503180d0         
                                g(  4) =         3.812970488946382d0    
                                  g(  5) =      3555.889503794543998d0  
             g(  6) =       284.029023525477442d0                       
                             g(  7) =       278.159456959177930d0       
                                  g(  8) =         0.467000314148584d0  
                          g(  9) =         0.755952269574951d0          
        g( 10) =         0.103963081719344d0                            
                  g( 11) =         0.000000000000000d0                  
          g( 12) =         0.000000000000000d0                          
                         g( 13) =         0.034604619091198d0           
           g( 14) =         0.447568930956009d0                         
                          g( 15) =       112.349387854963524d0          
                    g( 16) =         7.661018416622316d0                
       g( 17) =         0.000000000000000d0                             
          g( 18) =         7.000000000000000d0                          
                         g( 19) =        37.000000000000000d0           
                              g( 20) =         2.000000000000000d0      
               g( 21) =       354.000000000000000d0                     
                    g( 22) =         3.000000000000000d0                
                       g( 23) =         8.000000000000000d0             
                 g( 24) =         9.000000000000000d0                   
          g( 25) =         2.000000000000000d0                          
        g( 26) =         8.000000000000000d0                            
                      g( 27) =      1925.000000000000000d0              
                           g( 28) =        63.000000000000000d0         
       g( 29) =       281.000000000000000d0                             
                                    g( 30) =         6.000000000000000d0
                            g( 31) =        44.000000000000000d0        
                    g( 32) =     17294.000000000000000d0                
                                  g( 33) =        23.000000000000000d0  
               endif                                                    
                  if( i67 .eq. 3 )then                               
                    g(    1) =   1.00449566d0                           
                                      g(    2) =   0.01859477d0         
                       g(    3) =   0.97719994d0                        
                g(    4) =   1.07101045d0                               
                       g(    5) =   1345.89422853d0                     
          g(    6) =   185.25796624d0                                   
                                     g(    7) =   454.99070924d0        
                    g(    8) =   0.4030333d0                            
                                            g(    9) =   0.98286025d0   
        g(   10) =    0.0718894d0                                       
        g(   11) =    3.131d-05                                         
       g(   12) =    0.0d0                                              
           g(   13) =    0.20107243d0                                   
             g(   14) =    0.32253905d0                                 
                  g(   15) =    28.95573242d0                           
                               g(   16) =    4.35016847d0               
           g(   17) =    1.0d0                                          
                                    g(   18) =    59.94655842d0         
                         g(   19) =    10.40235209d0                    
         g(   20) =    45.5960256d0                                     
                   g(   21) =    5.33303639d0                           
                                          g(   22) =    11.0260158d0    
                  g(   23) =    5.62258459d0                            
                g(   24) =    11.8758362d0                              
                        g(   25) =    9.03405164d0                      
       g(   26) =    3.54369856d0                                       
                                 g(   27) =    1555.70240325d0          
              g(   28) =    3.27561139d0                                
                                g(   29) =    280.72701458d0            
                                g(   30) =    1.6899118d0               
               g(   31) =    7.72937078d0                               
                g(   32) =    2162.86269813d0                           
           g(   33) =    5.22324273d0                                   
                                            endif                       
       end                                                              
                        subroutine o37(k,i4,i32,w)                   
            implicit none                                               
                                 integer k,i32,j,i57,w                
                double precision i4                                     
                                       dimension i4(i32)                
                                   i57 = 0                            
              do j=1,k                                                  
         i57 = i57 + j                                              
                                                 enddo                  
                                                           do j=1,k     
                  i4(w+j-1) = dble(k-j+1)/dble(i57)                   
                                                   enddo                
         end                                                            
             subroutine o28(m,i8,k,g,i5,i2,i4,i32,i19,i49,i9,k3) 
                                            implicit none               
              integer m,i8,k,i32,i19,i49,i9,i,j,q,r,s,t,w             
       double precision i4(*),g(*),i5(*),i2(*),o25,i34,o16,k3(*),u,v
                              double precision z,k32,k33,i5i,i2i    
                                             u=i4(i9+1)                 
       v=i4(i9+2)                                                       
                       q=k-1                                            
      r=i49-1                                                         
                        s=i19-1                                         
                                          t=s+m                         
                                                           w = s-m      
                                k32 = k3(2)                         
                                          k33 = k3(3)               
                                                      do i=1,m          
                                              i34 = o25(i4)            
        z = i4(r+i) * o16(i4(1),i4(2))                                 
                       if( i34 .le. u )then                            
                z = z + i4(s+i)                                         
                                            else                        
              if( i34 .le. v )then                                     
                                          z = z + i4(t+i)               
                                     else                               
             do j=3,q                                                   
        if(i34.le.i4(i9+j)) goto 1                                     
                                                enddo                   
    1                                         z = z + i4(w+j*m+i)       
                         endif                                          
                                         endif                          
                                               i5i = i5(i)              
                         i2i = i2(i)                                    
                         if(z.lt.i5i)then                               
                           if(i34.ge.k32)then                        
              z = i5i + (i5i-z) * k33                                 
         if(z.gt.i2i) z = i2i                                           
                                                   else                 
                                                                 z = i5i
                                     endif                              
                                                                goto 2  
                 endif                                                  
                         if(z.gt.i2i)then                               
                             if(i34.ge.k32)then                      
         z = i2i - (z-i2i) * k33                                      
                                     if(z.lt.i5i) z = i5i               
                       else                                             
                                                  z = i2i               
                                                endif                   
                                endif                                   
    2                                   continue                        
                                               g(i) = z                 
                                 enddo                                  
                            if(i8.le.0) return                          
                                      do i=m-i8+1,m                     
                               g(i) = dnint(g(i))                       
                        enddo                                           
             if(i8.lt.m) return                                         
                                                           r = i19-1-m  
                              do j=1,k                                  
                    s = r+j*m                                           
                                         do i=1,m                       
                     t = s+i                                            
                                    if( g(i) .lt. i4(t) ) goto 88       
                               if( g(i) .gt. i4(t) ) goto 88            
                                   enddo                                
                        call o10(m,i8,g,i5,i2,i4,i32)        
                                                        return          
   88                                          continue                 
                                                               enddo    
                  end                                                   
                  subroutine k22(f,m,n,i0,l,x,g,i16)         
                                           implicit none                
              integer f,m,n,i0                                          
          double precision l(*),x(*),g(*),i16                           
      integer c,i,k23                                                 
                 double precision i17,k24,i56                       
                           if( f .le. 1 ) return                        
                                           if( n .le. 0 )then           
                                   i56 = l(1)                         
                                  k23 = 1                             
                                     do c=2,f                           
                              if( l(c) .lt. i56 )then                 
                        i56 = l(c)                                    
                k23 = c                                               
            endif                                                       
                                           enddo                        
                                                           endif        
                    if( n .ge. 1 )then                                  
                        call o31(i17,x(1),n,i0,i16)                    
                                        i56 = l(1)                    
                           k24 = i17                                  
                 k23 = 1                                              
                 do c=2,f                                               
                            call o31(i17,x((c-1)*n+1),n,i0,i16)        
           if( k24 .le. 0.0d0 )then                                   
                      if( i17 .le. 0.0d0 .and. l(c) .lt. i56 )then    
                                         i56 = l(c)                   
                      k24 = i17                                       
                               k23 = c                                
                             endif                                      
            else                                                        
                                             if( i17 .lt. k24 )then   
                            i56 = l(c)                                
              k24 = i17                                               
                                 k23 = c                              
        endif                                                           
                                               endif                    
                  enddo                                                 
                                            endif                       
                                  do c=1,f                              
       l(c) = l(k23)                                                  
             do i=1,n                                                   
                           x((c-1)*n+i) = x((k23-1)*n+i)              
                                                   enddo                
                do i=1,m                                                
               g((c-1)*m+i) = g((k23-1)*m+i)                          
                enddo                                                   
                                          enddo                         
                                                                    end 
          subroutine o9( g004, g019, io23,      
     &                           o, g005, io10,         
     &         io16)                                             
          implicit none                                                 
                            double precision g004, g019, io23
           double precision rio5,io16               
       integer o, g005,io10                             
                    io10 = 1                                     
                               io16 = 0.0d0                      
                         if( dabs( g019 - io23 ) .lt. 1.0d+32)then    
                if( dabs(g019-io23) .gt. 1.0d0 )then                  
                            if( o .le. 2 )then                          
             g004 = g004 - dsqrt(dabs(g019-io23))   
     &                   / (dble(g005))**0.1                   
                            else                                        
         g004 = g004 - dsqrt(dabs(g019-io23))       
             endif                                                      
                                                   else                 
        g004 = g004 - dabs(g019-io23)               
                                               endif                    
                io16 = rio5( g019, io23 )         
                                                                    else
                                g004 = g004 - 1.0d0   
             endif                                                      
                                                     end                
            subroutine o6( o,g008,io9 )  
                        implicit none                                   
                          integer i,o, i301,maxobji67               
                      double precision g008(*), io9        
          double precision g018, i17t, sung,io27            
          double precision fx25, fx26                
              do i=1,o                                                  
                                    g008(i) = 1.0d0                 
                                                                  enddo 
               if( i301(io9,0.0d0) .eq. 1 ) goto 999          
                             if( io9 .ge. 1.0d0 ) goto 999     
                        do i=1,o                                        
                       g008(i) = 0.0d0                              
                                      enddo                             
               maxobji67 = min(o,8)                                  
                                            do i=1,maxobji67         
               g018 = io9 * 10.0d0**dble(i-1)        
         g018 = g018 - io27(g018)       
               g018 = g018 * 10.0d0                 
          i17t = g018 - io27(g018)                
        if( i17t .gt. 0.99999999d0 ) g018 = g018 + 1
                 g018 = g018 - i17t                 
           g008(i) = g018                                 
                             enddo                                      
                                                       sung = 0.0d0     
                          do i=1,maxobji67                           
                       sung = sung + g008(i)                        
                                                    enddo               
                                  do i=1,maxobji67                   
                       g008(i) = g008(i) / sung                 
                                         enddo                          
                                    fx25 = 0.0d0               
                            do i=1,maxobji67                         
       if( g008(i) .gt. fx25 ) fx25 = g008(i) 
                                          enddo                         
                    fx26 = 1.0d0 / fx25              
               if( fx26 .gt. 1.0d0 )then                      
           do i=1,maxobji67                                          
                            g008(i) = g008(i) * fx26  
                        enddo                                           
        endif                                                           
  999               continue                                            
                                          end
                       subroutine o10(m,i8,g,i5,i2,i4,i32)   
                         implicit none                                  
                               integer m,i8,i32                         
                         double precision g(*),i5(*),i2(*),i4(*)        
                                  integer i,g017                     
                      double precision io6,o25                
                       io6 = 1.0d0 / dble(m)                  
                                     g017 = 0                        
                                       do i=1,m                         
                   if( i5(i) .ge. i2(i) ) goto 1                        
            if( o25(i4) .le. io6 )then                        
          if( g(i) .le. i5(i) )then                                     
                                    g(i) = i5(i) + 1.0d0                
                                                           g017 = 1  
                                 goto 1                                 
                                                      endif             
              if( g(i) .ge. i2(i) )then                                 
            g(i) = i2(i) - 1.0d0                                        
            g017 = 1                                                 
                                goto 1                                  
                                               endif                    
              if( o25(i4) .le. 0.5d0 )then                              
       g(i) = g(i) - 1.0d0                                              
                                                       g017 = 1      
             else                                                       
                                         g(i) = g(i) + 1.0d0            
                                               g017 = 1              
                         endif                                          
                                      endif                             
    1                                   continue                        
               enddo                                                    
                    if( g017 .eq. 0 ) call o30(m,i8,g,i5,i2,i4,i32)
                                                    end                 
                      subroutine io17( o,n,m,            
     & g014,                                              
     &   l,x,g,                                                         
     &                                      i17,i16,                    
     &                                g001,           
     &        io1,                                      
     &             g002,                               
     &        pl,k16,i306,i307,                                    
     &                            g008, io9,               
     &     io4, io20,                                  
     &                         io2,io10,           
     &                 io16, g005 )                     
                                                          implicit none 
            integer o, i                                                
      double precision g014,g001,i17,i16
         double precision l(*),x(*),g(*),g002(*)       
              double precision g008(*),pl(*),io1    
       double precision io9,k4t,fx19,i306(*),i307(*)    
       integer i449,i431,k,io4,best_i                 
      double precision i437,t1,t2,t3  
       double precision i305,fx07,t4,k6 
             double precision rio5,io16             
                  integer i432,i433,k16,n,m,i10ker, io20     
            integer io2,io10,g005         
              data i10ker /0/                                           
            if( io9 .ge. 1.0d0 ) goto 999                      
                                    if( i17 .gt. i16 ) goto 999         
         if( io1 .gt. i16 ) goto 999                    
                                               i449 = int(pl(1))       
                                    if( i449 .lt. 1 ) goto 999         
                                k6 = 1.0d-3                            
                                if( io4 .le. 1 )then      
          i10ker = 0                                                    
          io4 = 666                                       
                                                                  endif 
            i10ker = i10ker + 1                                         
                                                       i431 = 1+1    
        i432 = 1+k16*o+1                                           
                                     i433 = 1+k16*o+k16*n+1      
             if( io10 .eq. 1 .and. io16 .gt. 0.05d0 )then 
                               goto 888                                 
                                                       endif            
         i437 = 0.0d0                                                 
      t1  = 0.0d0                                                    
                                                        do i=1,o        
                                     if( g008(i) .gt. 0.0d0 )then   
                 if( l(i) .lt. g002(i) )then           
                                                       i437 = i437  
     &     + rio5( l(i), g002(i) )        
                                     endif                              
                        if( l(i) .gt. g002(i) )then    
                           t1  = t1                               
     &      + rio5( l(i), g002(i) )       
                                                 endif                  
                                         endif                          
                                  enddo                                 
             if( i437 .gt. 0.0d0 .and. t1 .le. 0.0d0 )then         
           if( g014 .ge. g001 )then     
               g014 = g001              
     &         - dabs(g001)*1.0d-8                    
                  endif                                                 
                                                          goto 888      
                               endif                                    
          if( i437 .le. 0.0d0 .and. t1 .gt. 0.0d0 )then            
            if( g014 .le. g001 )then    
           g014 = g001                  
     &                                          + t1                 
     &           * max( 1.0d0, dabs(g001))            
             endif                                                      
                        goto 999                                        
                                                             endif      
           if( i437 .le. 0.0d0 .and. t1 .le. 0.0d0 )then           
         if( g014 .gt. g001 )then       
                 g014 = g001            
          endif                                                         
                                                    goto 999            
                           endif                                        
                  if( io2 .eq. 0 .and. i449 .ge. 1 )then 
       if( g014 .le. g001 )then         
        g014 = g001                     
     &  + dabs(g001)*1.0d-8                           
                  goto 999                                              
                           endif                                        
                                           endif                        
        if( io2 .eq. 1 .and. i449 .eq. 1 )then           
              if( g014 .ge. g001 )then  
                    g014 = g001         
     &         - dabs(g001)*1.0d-8                    
                                 goto 999                               
                              endif                                     
                 endif                                                  
        if( io20 .ge. 2 .and.                                        
     &      g014 .gt. g001) goto 999    
                               if( io2 .eq. 1 ) goto 1    
         if( g014 .le. g001 )then       
                 g014 = g001            
     &                                   + k6                          
     &       * max( 1.0d0, dabs(g001))                
        endif                                                           
                                                          goto 999      
    1              continue                                             
                                          fx19=i305()     
                                                           do i=1,i449 
                                                 k4t = 0.0d0           
                            do k=1,o                                    
             k4t = k4t + rio5( g002(k), 
     &  pl(i431-1+ o*(i-1)+k) )                                      
        enddo                                                           
                                 if( k4t .lt. fx19 ) fx19 = k4t 
        if( k4t .le. k6 ) goto 2                                      
                                                      enddo             
              if( g014 .ge. g001 )then  
              g014 = g001               
     &                    - dabs(g001)*1.0d-8         
                                                          endif         
                                          goto 888                      
    2                                                   continue        
        call o1( o,                              
     &                l,                                                
     &       i17,i16, i306, i307,                                    
     &                          io9,                           
     &                                              0.0d0, g008,    
     &                            t2,                     
     &                                    g005)                
                                  call o1( o,    
     &                                      g002,      
     &       i17,i16, i306, i307,                                    
     &                                  io9,                   
     &                                 0.0d0, g008,                 
     &         t3,                                       
     &            g005)                                        
       if( g014 .gt. g001               
     &        .and. t2 .lt. t3 )then       
                         g014 = g001    
     &   - dabs(g001)*1.0d-8                          
               goto 888                                                 
           endif                                                        
            if( g014 .lt. g001          
     &         .and. t2 .ge. t3 )then      
       g014 = g001                      
     &                 + dabs(t2-t3)       
                           goto 999                                     
                                                          endif         
                                             if( i10ker .ge. 10 )then   
                                                       i10ker = 0       
                            goto 888                                    
                                                              endif     
  999           continue                                                
                 return                                                 
  888                continue                                           
                                 if(i449.lt.3) goto 999                
        fx07 = i305()                                  
                                            best_i = 1                  
                          do i=1,i449                                  
                                 call o1( o,     
     &                                     pl(i431-1+o*(i-1)+1),     
     &                                  i17,i16, i306, i307,         
     &           io9,                                          
     &                             0.0d0, g008,                     
     &                             t4,                    
     &      g005)                                              
      if( t4 .le. fx07 )then                        
                                fx07 = t4           
                                               best_i      = i          
                   endif                                                
                                enddo                                   
                              k4t = 0.0d0                              
            i = best_i                                                  
                                                               do k=1,o 
          k4t = k4t + rio5( g002(k),    
     &              pl(i431-1+ o*(i-1)+k) )                          
                                           enddo                        
                          if( k4t .le. 0.05d0 )then                    
                 else                                                   
                      do k=1,o                                          
                                                                enddo   
                              do k=1,o                                  
                 l(k) = pl(i431-1+ o*(best_i-1)+k)                   
                                                   enddo                
                                    do k=1,n                            
                x(k) = pl(i432-1+ n*(best_i-1)+k)                    
                        enddo                                           
              do k=1,m                                                  
                             g(k) = pl(i433-1+ m*(best_i-1)+k)       
                    enddo                                               
                  g014 = g001 - 1.0d-8  
     &                               * dabs(g001)     
                                                endif                   
               end                                                      
                 subroutine t5( g016, i67, p, z )                
                   implicit none                                        
                                     integer g016,i67,p,i,maxp       
                double precision z(*),io14(33),i5(33),i2(33)          
                      do i=1,33                                         
                                    io14(i)=0.0d0                     
                            enddo                                       
                 if( g016 .eq. 1 ) call io13( i67, io14 ) 
                 if( g016 .eq. 2 ) call io12( i67, io14 ) 
                                                   maxp = min(p,100)    
                                  do i=1,33                             
      z(i) = z(i) + io14(i) * (z(i)-z(i)/dsqrt(dble(maxp)))           
                                enddo                                   
                                                     i5(  1 ) = 0.1d0   
                          i2(  1 ) = 90.0d0                             
                               i5(  2 ) = 0.0d0                         
                                                i2(  2 ) = 1.0d0        
                              i5(  3 ) = 0.0d0                          
                                 i2(  3 ) = 1.0d0                       
                                                 i5(  4 ) = 0.0d0       
                                    i2(  4 ) = 9.0d0                    
                                          i5(  5 ) = 0.01d0             
                                                  i2(  5 ) = 9000.0d0   
      i5(  6 ) = 0.01d0                                                 
                          i2(  6 ) = 500.0d0                            
                        i5(  7 ) = 1.0d0                                
                                           i2(  7 ) = 5000.0d0          
                                 i5(  8 ) = 0.0d0                       
                          i2(  8 ) = 1.0d0                              
                            i5(  9 ) = 0.0d0                            
                         i2(  9 ) = 1.0d0                               
                 i5( 10 ) = 1.0d-12                                     
         i2( 10 ) = 0.5d0                                               
                                         i5( 11 ) = 0.0d0               
                                                      i2( 11 ) = 0.001d0
                              i5( 12 ) = 0.0d0                          
                                                   i2( 12 ) = 0.0d0     
                         i5( 13 ) = 0.0d0                               
                             i2( 13 ) = 0.999d0                         
                                           i5( 14 ) = 0.0d0             
                                   i2( 14 ) = 1.0d0                     
                                i5( 15 ) = 1.0d0                        
                                  i2( 15 ) = 500.0d0                    
                           i5( 16 ) = 4.0d0                             
                 i2( 16 ) = 12.0d0                                      
                 i5( 17 ) = 0.0d0                                       
                                                 i2( 17 ) = 0.0d0       
                      i5( 18 ) = 2.0d0                                  
            i2( 18 ) = 200.0d0                                          
                     i5( 19 ) = 1.0d0                                   
                                        i2( 19 ) = 1000.0d0             
              i5( 20 ) = 1.0d0                                          
                                                 i2( 20 ) = 100.0d0     
                                             i5( 21 ) = 1.0d0           
                                 i2( 21 ) = 500.0d0                     
                                                       i5( 22 ) = 1.0d0 
        i2( 22 ) = 12.0d0                                               
                     i5( 23 ) = 1.0d0                                   
                      i2( 23 ) = 12.0d0                                 
                                       i5( 24 ) = 1.0d0                 
                    i2( 24 ) = 12.0d0                                   
                                                       i5( 25 ) = 1.0d0 
                                                   i2( 25 ) = 12.0d0    
                                                  i5( 26 ) = 1.0d0      
                              i2( 26 ) = 500.0d0                        
      i5( 27 ) = 2.0d0                                                  
                        i2( 27 ) = 2000.0d0                             
                                                  i5( 28 ) = 2.0d0      
                    i2( 28 ) = 100.0d0                                  
                    i5( 29 ) = 2.0d0                                    
                            i2( 29 ) = 2000.0d0                         
       i5( 30 ) = 1.0d0                                                 
                                            i2( 30 ) = 12.0d0           
                   i5( 31 ) = 1.0d0                                     
                                           i2( 31 ) = 100.0d0           
                                             i5( 32 ) = 10.0d0          
                               i2( 32 ) = 90000.0d0                     
                                                i5( 33 ) = 1.0d0        
                                                  i2( 33 ) = 30.0d0     
      do i=1,33                                                         
                        if(z(i).lt.i5(i)) z(i)=i5(i)                    
                    if(z(i).gt.i2(i)) z(i)=i2(i)                        
                                                          enddo         
                                                   end                  
                           subroutine io13(i67,g)           
                                                         implicit none  
         integer i67                                                 
                         double precision g(*)                          
                         if( i67 .eq . 1 )then                       
                                               g(  1) =  4.29154553d0   
                                          g(  2) =  1.47194224d0        
                                               g(  3) =  1.54997898d0   
                        g(  4) =  -0.17704074d0                         
                                                  g(  5) =  -1.0d0      
               g(  6) =  0.56172152d0                                   
                g(  7) =  0.05515924d0                                  
                    g(  8) =  -0.3771567d0                              
                                       g(  9) =  3.74633499d0           
                  g( 10) =  1.78936255d0                                
        g( 11) =  -0.15518801d0                                         
         g( 12) =  -0.98360507d0                                        
                        g( 13) =  0.58322553d0                          
                              g( 14) =  -0.9608422d0                    
       g( 15) =  1.02212226d0                                           
                                     g( 16) =  3.58965909d0             
                  g( 17) =  2.13785201d0                                
                                      g( 18) =  0.10303835d0            
                                     g( 19) =  4.04614941d0             
                         g( 20) =  0.05596182d0                         
                                g( 21) =  -0.74156199d0                 
                                       g( 22) =  5.19929225d0           
                            g( 23) =  4.10382737d0                      
                                g( 24) =  4.4548273d0                   
                              g( 25) =  4.41742088d0                    
                      g( 26) =  0.38736892d0                            
               g( 27) =  0.17670075d0                                   
           g( 28) =  4.15725493d0                                       
                                   g( 29) =  -0.22783536d0              
                                                  g( 30) =  0.58237005d0
                             g( 31) =  3.04829469d0                     
                                                    g( 32) =  -1.0d0    
            g( 33) =  -0.28213036d0                                     
                        endif                                           
                     if( i67 .eq . 2 )then                           
                                        g(  1) =   1.07767821d0         
                         g(  2) =   -0.44761809d0                       
                    g(  3) =   0.65368561d0                             
                                   g(  4) =   -1.0d0                    
                      g(  5) =   -0.73027076d0                          
      g(  6) =   -0.96722505d0                                          
                                              g(  7) =   -0.16648158d0  
                     g(  8) =   -0.9630052d0                            
       g(  9) =   -0.43043102d0                                         
                   g( 10) =    -0.33013524d0                            
                                   g( 11) =    -0.53458138d0            
        g( 12) =    0.04519436d0                                        
                                        g( 13) =    0.21059818d0        
                             g( 14) =    1.85997072d0                   
           g( 15) =    1.00887824d0                                     
               g( 16) =    0.05296895d0                                 
                                          g( 17) =    -0.75654657d0     
               g( 18) =    3.99405505d0                                 
        g( 19) =    2.10079382d0                                        
            g( 20) =    0.71820448d0                                    
                       g( 21) =    -0.79308107d0                        
                                g( 22) =    0.71794093d0                
              g( 23) =    1.59649466d0                                  
                                             g( 24) =    -0.42941671d0  
                                     g( 25) =    0.05688063d0           
                              g( 26) =    -0.36715641d0                 
                                                g( 27) =    0.10866524d0
        g( 28) =    -0.40182486d0                                       
                     g( 29) =    0.87367028d0                           
                       g( 30) =    -0.15730903d0                        
                                         g( 31) =    2.53653413d0       
                                g( 32) =    -0.2767428d0                
                              g( 33) =    -0.21904374d0                 
                                                     endif              
                                                   end                  
                               subroutine io12(i67,g)       
                                                           implicit none
                 integer i67                                         
                                       double precision g(*)            
                                             if( i67 .eq . 1 )then   
                                g(  1) =  2.08181871d0                  
                                                 g(  2) =  0.19941008d0 
                 g(  3) =  0.663559d0                                   
                              g(  4) =  0.2297929d0                     
                                      g(  5) =  0.45194597d0            
               g(  6) =  -0.73231382d0                                  
              g(  7) =  1.12731735d0                                    
          g(  8) =  0.15438953d0                                        
                           g(  9) =  0.92589979d0                       
                                      g( 10) =   0.66252956d0           
                    g( 11) =   1.51928446d0                             
                 g( 12) =   -0.42972569d0                               
                   g( 13) =   1.67373601d0                              
                                            g( 14) =   1.88837352d0     
                        g( 15) =   -0.21119113d0                        
                     g( 16) =   0.02374419d0                            
                                  g( 17) =   0.45582854d0               
                         g( 18) =   -0.27602405d0                       
             g( 19) =   -0.53242829d0                                   
                  g( 20) =   0.14369251d0                               
                                           g( 21) =   -0.00105148d0     
                   g( 22) =   -0.77738419d0                             
              g( 23) =   -0.53794328d0                                  
                      g( 24) =   -0.55895883d0                          
                                             g( 25) =   0.28001733d0    
                          g( 26) =   -0.11034496d0                      
                           g( 27) =   -0.40392603d0                     
                                            g( 28) =   0.29138653d0     
                                            g( 29) =   0.2377215d0      
                                                g( 30) =   0.20211495d0 
                          g( 31) =   0.138778d0                         
       g( 32) =   0.65156465d0                                          
                         g( 33) =   -0.99249717d0                       
                                                              endif     
                      if( i67 .eq . 2 )then                          
                         g(  1) =   -0.4312276d0                        
                                     g(  2) =   -0.15189541d0           
                              g(  3) =   0.03960097d0                   
                         g(  4) =   -0.33023865d0                       
                          g(  5) =   0.37549423d0                       
                           g(  6) =   0.46619011d0                      
                       g(  7) =   0.01955984d0                          
                     g(  8) =   -0.61703271d0                           
                g(  9) =   -0.10107408d0                                
                                                g( 10) =    0.0082648d0 
                             g( 11) =    -0.03994651d0                  
                          g( 12) =    -0.21620967d0                     
                                      g( 13) =    0.12164474d0          
                                        g( 14) =    0.12107576d0        
                                      g( 15) =    0.18700566d0          
       g( 16) =    -0.56889807d0                                        
                   g( 17) =    0.00153063d0                             
          g( 18) =    -0.07262216d0                                     
                         g( 19) =    0.80049358d0                       
                                         g( 20) =    -0.20426258d0      
                       g( 21) =    0.32111134d0                         
                        g( 22) =    0.20215425d0                        
                                    g( 23) =    0.25647943d0            
                             g( 24) =    0.1400869d0                    
                       g( 25) =    0.16738659d0                         
                      g( 26) =    1.03552112d0                          
                    g( 27) =    0.19846093d0                            
           g( 28) =    -0.06430121d0                                    
                                        g( 29) =    0.12672987d0        
                g( 30) =    -0.33030999d0                               
                                            g( 31) =    0.15749225d0    
                                              g( 32) =    -0.05867185d0 
                                          g( 33) =    0.59177165d0      
                                 endif                                  
                                                                     end
                 integer function i301( a, b )                         
                              implicit none                             
                                       double precision a, b            
                                  i301 = 0                             
                         if( a.lt.b ) return                            
                                    if( a.gt.b ) return                 
                                     i301 = 1                          
                         end                                            
              integer function i304(g)                                
         implicit none                                                  
                                 double precision g                     
                                                     i304 = 0         
              if( g .ne. g ) i304 = 1                                 
                          end                                           
              double precision function i305()               
                                                  implicit none         
                   double precision k26,i44,k25             
                      integer k27,k28,i                        
                                        data k27 /1/              
                                 data k26 /1.0d0/                  
                       if( k27 .eq. 1)then                        
                             k27    = 0                           
                               i44  = 0.0d0                           
                           k28 = 0                                   
                        do i = 1,99                                     
               i44 = k26                                         
                                        k26 = k26 * 10.0d0    
                        k28 = k28 + 1                             
            if( i44 .ge. k26 ) goto 1                            
                                  enddo                                 
    1  continue                                                         
                                 k25 = dble(k28-1)              
          k26 = dabs( 10.0d0**k25 )                           
                 if( k26 .lt. 1.0d+16) k26 = 1.0d+16          
                          endif                                         
                                  i305 = k26            
                                          end
               subroutine io7(l,x,o,n)                   
                       implicit none                                    
                integer o,n,i,i304                                    
                 double precision l(*),x(*), i305, k26  
                      k26 = i305()                      
                   do i=1,o                                             
                   if( i304(l(i)) .eq. 1 ) l(i) = k26            
          if( l(i) .gt. k26 )  l(i) = k26                     
                                            enddo                       
                        do i=1,n                                        
                        if( i304(x(i)) .eq. 1 )  x(i) = - k26    
                      if( x(i) .lt. - k26 ) x(i) = - k26      
                                          enddo                         
                                                             end        
                              subroutine io7_lomfy(l,o)  
             implicit none                                              
                                      integer o,i,i304                
                   double precision l(*), i305, k26     
              k26 = i305()                              
                        do i=1,o                                        
            if( i304(l(i)) .eq. 1 ) l(i) = k26                   
                          if( l(i) .gt. k26 )  l(i) = k26     
                            enddo                                       
             end                                                        
          double precision function rio5(a,b)              
                                            implicit none               
         double precision a,b,maxab                                     
                                 maxab = max( dabs(a), dabs(b) )        
                          rio5 = dabs(a-b)                 
       rio5 = rio5 / max( 1.0d0, maxab )      
               end                                                      
                               subroutine o18(i67s,i15)      
         implicit none                                                  
                                          character i15*60              
                                integer i67s(*),i                    
         do i = 1,60                                                    
        call alphabet(i15(i:i),i67s(i))                              
                                                enddo                   
                                                         end            
c end                  ol code                                          
      subroutine alphabet(a,b)      
      implicit none      
      character a*1
      integer b      
      b = 0      
      if(a(1:1).eq.'A') b = 52
      if(a(1:1).eq.'B') b = 28
      if(a(1:1).eq.'C') b = 49
      if(a(1:1).eq.'D') b = 30
      if(a(1:1).eq.'E') b = 31
      if(a(1:1).eq.'F') b = 32
      if(a(1:1).eq.'G') b = 33
      if(a(1:1).eq.'H') b = 34
      if(a(1:1).eq.'I') b = 35
      if(a(1:1).eq.'J') b = 36
      if(a(1:1).eq.'K') b = 37
      if(a(1:1).eq.'L') b = 38
      if(a(1:1).eq.'M') b = 39
      if(a(1:1).eq.'N') b = 40
      if(a(1:1).eq.'O') b = 41
      if(a(1:1).eq.'P') b = 42      
      if(a(1:1).eq.'Q') b = 43
      if(a(1:1).eq.'R') b = 44
      if(a(1:1).eq.'S') b = 45
      if(a(1:1).eq.'T') b = 46
      if(a(1:1).eq.'U') b = 47
      if(a(1:1).eq.'V') b = 48
      if(a(1:1).eq.'W') b = 29
      if(a(1:1).eq.'X') b = 50
      if(a(1:1).eq.'Y') b = 51
      if(a(1:1).eq.'Z') b = 27      
      if(a(1:1).eq.'0') b = 53
      if(a(1:1).eq.'1') b = 54
      if(a(1:1).eq.'2') b = 55
      if(a(1:1).eq.'3') b = 56
      if(a(1:1).eq.'4') b = 57
      if(a(1:1).eq.'5') b = 58
      if(a(1:1).eq.'6') b = 59
      if(a(1:1).eq.'7') b = 60
      if(a(1:1).eq.'8') b = 61  
      if(a(1:1).eq.'9') b = 62              
      if(a(1:1).eq.'a') b = 23
      if(a(1:1).eq.'b') b = 2
      if(a(1:1).eq.'c') b = 3
      if(a(1:1).eq.'d') b = 16
      if(a(1:1).eq.'e') b = 5
      if(a(1:1).eq.'f') b = 13
      if(a(1:1).eq.'g') b = 7
      if(a(1:1).eq.'h') b = 8
      if(a(1:1).eq.'i') b = 9
      if(a(1:1).eq.'j') b = 10
      if(a(1:1).eq.'k') b = 11
      if(a(1:1).eq.'l') b = 12
      if(a(1:1).eq.'m') b = 6
      if(a(1:1).eq.'n') b = 14
      if(a(1:1).eq.'o') b = 15
      if(a(1:1).eq.'p') b = 4
      if(a(1:1).eq.'q') b = 17
      if(a(1:1).eq.'r') b = 18
      if(a(1:1).eq.'s') b = 19
      if(a(1:1).eq.'t') b = 20
      if(a(1:1).eq.'u') b = 21
      if(a(1:1).eq.'v') b = 22
      if(a(1:1).eq.'w') b = 1
      if(a(1:1).eq.'x') b = 24
      if(a(1:1).eq.'y') b = 25
      if(a(1:1).eq.'z') b = 26                   
      if(a(1:1).eq.'_') b = 64  
      if(a(1:1).eq.'(') b = 65
      if(a(1:1).eq.')') b = 66
      if(a(1:1).eq.'+') b = 67
      if(a(1:1).eq.'-') b = 68
      if(a(1:1).eq.'&') b = 69 
      if(a(1:1).eq.'.') b = 70
      if(a(1:1).eq.',') b = 71
      if(a(1:1).eq.':') b = 72
      if(a(1:1).eq.';') b = 73 
      if(a(1:1).eq.'*') b = 74
      if(a(1:1).eq.'=') b = 75
      if(a(1:1).eq.'/') b = 76                           
      if(a(1:1).eq.'!') b = 80  
      if(a(1:1).eq.'[') b = 83
      if(a(1:1).eq.']') b = 84            
      end
       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine handles all printing commands for MIDACO.
C     This subroutine will also check the MAXEVAL and MAXTIME criteria.
C     Note that this subroutine is called independently from MIDACO and
C     MIDACO itself does not include any print commands (due to 
C     compiler portability and robustness).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MIDACO_PRINT(C,PRINTEVAL,SAVE2FILE,IFLAG,ISTOP,F,G,X,
     &                        XL,XU,O,N,NI,M,ME,RW,PF,MAXEVAL,MAXTIME,
     &                        PARAM,P,T,KEY)           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Increase size, if problems with N or M+2*O > 10000 are solved    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION BESTX(10000), BESTG(10000)      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER C,PRINTEVAL,SAVE2FILE,IFLAG,ISTOP,EVAL,O,N,NI,M,ME,I,P,T
      INTEGER TIC,Q,IOUT1,IOUT2,IOUT3,IOUT4,MAXEVAL,MAXTIME,UPDATE,PSIZE
      INTEGER KF,KG,KRES,KX,WF,WG,WRES,WX,KBEST,WBEST,PFMAX,same

      integer EXTRAOFFSET

      DOUBLE PRECISION TNOW,TSTART,TMAX,RW(*),BESTF,BESTR,F(*),G(*),X(*)
      DOUBLE PRECISION XL(*),XU(*),ACC,PARAM(*),DUMMY_F,DUMMY_VIO,PF(*)
      CHARACTER KEY*60
      DATA KF,KG,KRES,KX,WF,WG,WRES,WX,PFMAX /0,0,0,0,0,0,0,0,0/      
      DATA KBEST,WBEST,TIC,Q,EVAL,IOUT1,IOUT2,IOUT3 /0,0,0,0,0,0,0,0/   
      DATA BESTF,BESTR,DUMMY_F,DUMMY_VIO /0.0D0,0.0D0,0.0D0,0.0D0/
      DATA TNOW,TSTART,TMAX,ACC /0.0D0,0.0D0,0.0D0,0.0D0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fortran printing units, freely change those values
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IOUT1 = 36 ! MIDACO_SCREEN.TXT
      IOUT2 = 37 ! MIDACO_SOLUTION.TXT
      IOUT3 = 38 ! MIDACO_PARETOFRONT.TXT
      IOUT4 = 39 ! MIDACO_HISTORY.TXT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(C.EQ.1)THEN             
          IFLAG = 0
          ISTOP = 0      
          TMAX = DBLE(MAXTIME)           
          CALL GET_TIME(TSTART,P,T)        
          EVAL = 0                
          IF(PARAM(1).LE.0.0D0)THEN
              ACC = 1.0D-3
          ELSE
              ACC = PARAM(1)
          ENDIF              
C         Note that this offset differs      
          EXTRAOFFSET = 5*(P+O+N+M)+100
          Q    = 102*N+(M+2*O)+516 + EXTRAOFFSET
          KX   = 9          
          KF   = 9+1+N
          KG   = 9+1+N
          KRES = 9+1+N+1+ M
          IF( O .GT. 1 ) KRES = 9+1+N+1+ (M+2*O)
          WX   = Q
          WF   = Q+1+N
          WG   = Q+1+N
          WRES = Q+1+N+1+ M
          IF( O .GT. 1 ) WRES = Q+1+N+1+ (M+2*O)            
          IF( SAVE2FILE.GE.1 .AND. PRINTEVAL.GE.1 )THEN          
              OPEN(IOUT1,FILE='MIDACO_SCREEN.TXT'  ,STATUS='UNKNOWN')
              OPEN(IOUT2,FILE='MIDACO_SOLUTION.TXT',STATUS='UNKNOWN') 
          ENDIF     
          BESTF = 1.0D+99
          BESTR = 1.0D+99 
          DUMMY_F   = 1.0D+99
          DUMMY_VIO = 1.0D+99      
          TIC = 0             
          PFMAX = 1000
          IF(PARAM(10).GE. 1.0D0) PFMAX =  NINT(PARAM(10))
          IF(PARAM(10).LE.-1.0D0) PFMAX = -NINT(PARAM(10))          
          IF(PRINTEVAL.GE.1)THEN
              CALL PRINT_HEAD( P, O, N,NI,M,ME, PARAM, MAXEVAL, MAXTIME,
     &                         PRINTEVAL, SAVE2FILE, KEY, 0)
              IF(SAVE2FILE.GE.1)THEN
                CALL PRINT_HEAD( P, O, N,NI,M,ME, PARAM,MAXEVAL,MAXTIME,
     &                        PRINTEVAL, SAVE2FILE, KEY, IOUT1)
              ENDIF            
              IF(SAVE2FILE.GE.1) WRITE(IOUT2,3)
    3         FORMAT('MIDACO - SOLUTION',/,
     &           '-----------------',/,
     &'This file saves the current best solution X found by MIDACO.',/,
     &'This file is updated after every PRINTEVAL function evaluation,',
     &/,'if X has been improved.',/,/)              
              IF(SAVE2FILE.GE.1) CALL FLUSH_OUTPUT( IOUT1 )
              IF(SAVE2FILE.GE.1) CALL FLUSH_OUTPUT( IOUT2 )   
          ENDIF                                         
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
      IF(SAVE2FILE.GE.2) CALL SAVE_HISTORY(SAVE2FILE,P,O,N,M,F,G,X,
     &                                     IOUT4,IFLAG,ISTOP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      IF(C.EQ.2)THEN  
          CALL GET_TIME(TNOW,P,T)
          TNOW = TNOW - TSTART 
          EVAL = EVAL + P                               
          IF(IFLAG.GE.10)THEN          
              CALL WARNINGS_AND_ERRORS( IFLAG, 0 )              
              IF(SAVE2FILE.GE.1)THEN
                  CALL WARNINGS_AND_ERRORS( IFLAG, IOUT1 )
                  CALL WARNINGS_AND_ERRORS( IFLAG, IOUT2 )           
              ENDIF  
          ENDIF  
          IF(PRINTEVAL.GE.1)THEN          
              TIC = TIC + P
              IF(TIC.GE.PRINTEVAL.OR.EVAL.EQ.P.OR.IFLAG.GE.1)THEN
                  IF(EVAL.LT.0)THEN
                    TIC = 0
                    EVAL = -EVAL-2*PRINTEVAL
                  ENDIF
                  IF(EVAL.GT.P) TIC = 0      
                  IF(same(RW(KRES),RW(WRES)).EQ.1)THEN            
                    KBEST=KF
                    WBEST=WF
                  ELSE                 
                    KBEST=KRES
                    WBEST=WRES
                  ENDIF                 
                  IF(RW(WBEST).LT.RW(KBEST).OR.
     &              (IFLAG.GE.1.OR.IFLAG.EQ.-300.OR.IFLAG.EQ.-500))THEN                  
                      BESTF = RW(WF)
                      BESTR = RW(WRES)
                      IF( O .LE. 1 )THEN                      
                        DO I=1,M
                          BESTG(I) = RW(WG+I)                      
                        ENDDO
                      ELSE
                        DO I=1,M+2*O
                          BESTG(I) = RW(WG+I)                     
                        ENDDO
                      ENDIF
                      DO I=1,N
                        BESTX(I) = RW(WX+I)
                      ENDDO
                  ELSE                  
                      BESTF = RW(KF)
                      BESTR = RW(KRES)
                      IF( O .LE. 1 )THEN
                        DO I=1,M
                          BESTG(I) = RW(KG+I)                    
                        ENDDO      
                      ELSE
                        DO I=1,M+2*O
                          BESTG(I) = RW(KG+I)                    
                        ENDDO
                      ENDIF
                      DO I=1,N
                        BESTX(I) = RW(KX+I)
                      ENDDO                      
                  ENDIF                               
                  PSIZE = NINT(PF(1))   
                  IF(IFLAG.LT.100)THEN   
                  CALL PRINT_LINE(O,EVAL,TNOW,BESTF,BESTR,PSIZE,0)    
                  IF(SAVE2FILE.GE.1)THEN                    
                    CALL PRINT_LINE(O,EVAL,TNOW,BESTF,BESTR,PSIZE,IOUT1)
                  ENDIF               
                  IF(SAVE2FILE.GE.1)THEN   
                    UPDATE = 0                 
                    IF( (BESTR.LT.DUMMY_VIO) .OR.
     &                  (same(BESTR,DUMMY_VIO).EQ.1
     &                   .AND.BESTF.LT.DUMMY_F) )THEN 
                       DUMMY_F   = BESTF
                       DUMMY_VIO = BESTR
                       UPDATE    = 1 
                    ENDIF                                   
                    IF(UPDATE.EQ.1)THEN
                     WRITE(IOUT2,31)
   31                FORMAT(/,/,'            CURRENT BEST SOLUTION') 
                     CALL PRINT_SOLUTION( O,N,M,ME, BESTX, BESTG, BESTF,
     &                    BESTR, PF, XL,XU,ACC,EVAL, TNOW, IFLAG, IOUT2)
                     CONTINUE     
                    ENDIF                    
                    CALL FLUSH_OUTPUT( IOUT1 )
                    CALL FLUSH_OUTPUT( IOUT2 )
                  ENDIF  
                  IF( O .GT. 1 )THEN
                    IF(SAVE2FILE.GE.1)  CALL PRINT_PARETOFRONT( O,M,N,
     &                                       PF,BESTX,BESTG,PFMAX,IOUT3)
                  ENDIF
                  ENDIF                        
              ENDIF  
          ENDIF    
          IF(ISTOP.EQ.0)THEN
              IF(TNOW.GE.TMAX     ) IFLAG = -999
              IF(EVAL.GE.MAXEVAL-1)THEN            
                IF(MAXEVAL.LT.99999999) IFLAG = -999
              ENDIF
C             Special Case: maxeval = 1         
              IF( MAXEVAL.LE.1 )THEN
                IFLAG=0
                ISTOP=1
                DO I=1,N
                 X(I) = BESTX(I)
                ENDDO
              ENDIF              
          ENDIF       
          IF(ISTOP.GE.1)THEN
            DO I=1,M
              BESTG(I)=G(I)
            ENDDO
            IF(O.GE.2)THEN
              BESTF = RW(WF)
              DO I=1,O*2
                BESTG(M+I)=RW(WG+M+I) 
              ENDDO
            ENDIF
            if(O.EQ.1) bestf = F(1)
            IF(PRINTEVAL.GT.0)THEN
            CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &               O,N,M,ME,X,BESTG,BESTF,XL,XU,RW,ACC,WRES,PF,PARAM, 
     &                        0 )     
            IF(SAVE2FILE.GT.0)THEN
               CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &               O,N,M,ME,X,BESTG,BESTF,XL,XU,RW,ACC,WRES,PF,PARAM, 
     &                           IOUT1 )
               CALL PRINT_FINAL( IFLAG,TNOW,TMAX,EVAL,MAXEVAL,
     &               O,N,M,ME,X,BESTG,BESTF,XL,XU,RW,ACC,WRES,PF,PARAM, 
     &                           IOUT2 )            
               CALL FLUSH_OUTPUT( IOUT1 )
               CALL FLUSH_OUTPUT( IOUT2 )
               CALL CLOSE_OUTPUT( IOUT1 )
               CALL CLOSE_OUTPUT( IOUT2 )              
            ENDIF
            ENDIF
          ENDIF
          RETURN  
      ENDIF  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRINT_HEAD( P, O, N,NI,M,ME, PAR, MAXEVAL, MAXTIME, 
     &                       PRINTEVAL, SAVE2FILE, KEY, IOUT)     
      IMPLICIT NONE
      INTEGER P,O,N,NI,M,ME,IOUT,MAXTIME,MAXEVAL,SAVE2FILE,PRINTEVAL,I
      DOUBLE PRECISION PAR(*)
      CHARACTER KEY*60 
      WRITE(IOUT,1) KEY,O,P,N,MAXEVAL,NI,
     &              MAXTIME,M,PRINTEVAL,ME,SAVE2FILE
      DO I = 1,13
          IF(PAR(I).GT.0.0D0 .OR. PAR(I).LT.0.0D0) GOTO 66
      ENDDO
      IF(PRINTEVAL.GE.1) WRITE(IOUT,77)     
      GOTO 67
   66 IF(PAR( 1).GT.0.0D0.OR.PAR( 1).LT.0.0D0) WRITE(IOUT,201) PAR( 1)
      IF(PAR( 2).GT.0.0D0.OR.PAR( 2).LT.0.0D0) WRITE(IOUT,202) PAR( 2)
      IF(PAR( 3).GT.0.0D0.OR.PAR( 3).LT.0.0D0) WRITE(IOUT,203) PAR( 3)
      IF(PAR( 4).GT.0.0D0.OR.PAR( 4).LT.0.0D0) WRITE(IOUT,204) PAR( 4)
      IF(PAR( 5).GT.0.0D0.OR.PAR( 5).LT.0.0D0) WRITE(IOUT,205) PAR( 5)
      IF(PAR( 6).GT.0.0D0.OR.PAR( 6).LT.0.0D0) WRITE(IOUT,206) PAR( 6)
      IF(PAR( 7).GT.0.0D0.OR.PAR( 7).LT.0.0D0) WRITE(IOUT,207) PAR( 7)
      IF(PAR( 8).GT.0.0D0.OR.PAR( 8).LT.0.0D0) WRITE(IOUT,208) PAR( 8)
      IF(PAR( 9).GT.0.0D0.OR.PAR( 9).LT.0.0D0) WRITE(IOUT,209) PAR( 9)
      IF(PAR(10).GT.0.0D0.OR.PAR(10).LT.0.0D0) WRITE(IOUT,210) PAR(10)
      IF(PAR(11).GT.0.0D0.OR.PAR(11).LT.0.0D0) WRITE(IOUT,211) PAR(11)
      IF(PAR(12).GT.0.0D0.OR.PAR(12).LT.0.0D0) WRITE(IOUT,212) PAR(12)
      IF(PAR(13).GT.0.0D0.OR.PAR(13).LT.0.0D0) WRITE(IOUT,213) PAR(13)                     
      WRITE(IOUT,188)          
   67 CONTINUE      
      IF( O .EQ. 1 ) WRITE(IOUT,302)
      IF( O .GT. 1 ) WRITE(IOUT,303)          
    1 FORMAT(/,
     &       ' MIDACO 6.0    (www.midaco-solver.com)',/,
     &       ' -------------------------------------',/,/,
     &       ' LICENSE-KEY:  ',A60,/,/,           
     &       ' ----------------------------------------',/,
     &       ' | OBJECTIVES',I5,' | PARALLEL',I10,' |',/,
     &       ' |--------------------------------------|',/,     
     &       ' | N', I14, ' | MAXEVAL',    I11,' |',/,
     &       ' | NI',I13, ' | MAXTIME',    I11,' |',/,     
     &       ' | M', I14, ' | PRINTEVAL',  I9,' |',/,     
     &       ' | ME',I13 ,' | SAVE2FILE',  I9,' |',/,
     &       ' |--------------------------------------|')  
   77 FORMAT(' | PARAMETER:    All by default (0)     |',/,
     &       ' ----------------------------------------')   
  201 FORMAT(' | PARAM( 1) ',E14.7,' ACCURACY    |')
  202 FORMAT(' | PARAM( 2) ',E14.7,' RANDOM-SEED |')    
  203 FORMAT(' | PARAM( 3) ',E14.7,' FSTOP       |') 
  204 FORMAT(' | PARAM( 4) ',E14.7,' ALGOSTOP    |')     
  205 FORMAT(' | PARAM( 5) ',E14.7,' EVALSTOP    |') 
  206 FORMAT(' | PARAM( 6) ',E14.7,' FOCUS       |')   
  207 FORMAT(' | PARAM( 7) ',E14.7,' ANTS        |')  
  208 FORMAT(' | PARAM( 8) ',E14.7,' KERNEL      |')
  209 FORMAT(' | PARAM( 9) ',E14.7,' ORACLE      |')
  210 FORMAT(' | PARAM(10) ',E14.7,' PARETOMAX   |') 
  211 FORMAT(' | PARAM(11) ',E14.7,' EPSILON     |')   
  212 FORMAT(' | PARAM(12) ',E14.7,' BALANCE     |')
  213 FORMAT(' | PARAM(13) ',E14.7,' CHARACTER   |')    
  188 FORMAT(' ----------------------------------------') 
  302   FORMAT(/,/,
     &  ' [     EVAL,    TIME]        OBJECTIVE FUNCTION VALUE    ',
     &  '     VIOLATION OF G(X)',/,
     &  ' -----------------------------------------------',
     &  '-------------------------------') 
  303   FORMAT(/,/,
     &  ' [     EVAL,    TIME]   MULTI-OBJECTIVE PROGRESS  ',     
     &  ' VIOLATION OF G(X)',/,
     &  ' --------------------------------------------------',
     &  '-----------------   [PARETO]')                 
      END
      SUBROUTINE PRINT_LINE( O, EVAL, TNOW, F, VIO, PSIZE, IOUT)
      IMPLICIT NONE
      INTEGER O, EVAL, PSIZE, IOUT
      DOUBLE PRECISION TNOW, F, VIO
      IF( O .EQ. 1 )THEN
        IF(DABS(F).LE.1.0D+10)THEN
          IF(VIO.LE.1.0D+5)THEN
             WRITE(IOUT,1) EVAL, TNOW, F, VIO           
          ELSE
             IF(O.EQ.1) WRITE(IOUT,2) EVAL, TNOW, F, VIO    
          ENDIF
        ELSE
          IF(VIO.LE.1.0D+5)THEN
             WRITE(IOUT,3) EVAL, TNOW, F, VIO          
          ELSE
             WRITE(IOUT,4) EVAL, TNOW, F, VIO       
          ENDIF            
        ENDIF             
      ELSE
        IF(DABS(F).LE.1.0D+10)THEN
          IF(VIO.LE.1.0D+5)THEN
             WRITE(IOUT,111) EVAL, TNOW, F, VIO, PSIZE     
          ELSE
             WRITE(IOUT,112) EVAL, TNOW, F, VIO, PSIZE      
          ENDIF
        ELSE
          IF(VIO.LE.1.0D+5)THEN
             WRITE(IOUT,113) EVAL, TNOW, F, VIO, PSIZE       
          ELSE
             WRITE(IOUT,114) EVAL, TNOW, F, VIO, PSIZE 
          ENDIF            
        ENDIF
      ENDIF           
    1 FORMAT(
     &' [',I9,',',F8.2,']        F(X):',F19.8,'         VIO:',F13.6) 
    2 FORMAT(
     &' [',I9,',',F8.2,']        F(X):',F19.8,'         VIO:',D13.6)
    3 FORMAT(
     &' [',I9,',',F8.2,']        F(X):',E19.10,'         VIO:',F13.6)
    4 FORMAT(
     &' [',I9,',',F8.2,']        F(X):',E19.10,'         VIO:',D13.6)   
  111 FORMAT(
     &' [',I9,',',F8.2,']   PRO:', F20.8,'   VIO:',F13.6,'   [',I6,']')
  112 FORMAT(
     &' [',I9,',',F8.2,']   PRO:', F20.8,'   VIO:',D13.6,'   [',I6,']')
  113 FORMAT(
     &' [',I9,',',F8.2,']   PRO:',E20.10,'   VIO:',F13.6,'   [',I6,']')
  114 FORMAT(
     &' [',I9,',',F8.2,']   PRO:',E20.10,'   VIO:',D13.6,'   [',I6,']') 
      END
      SUBROUTINE PRINT_SOLUTION( O, N, M, ME, X, G, F, VIO, PF,
     &                           XL, XU, ACC, EVAL, TNOW, IFLAG, IOUT)  
      IMPLICIT NONE
      INTEGER O, N,M,ME,EVAL,IFLAG,IOUT,I,J,same
      DOUBLE PRECISION X(*),G(*),F,VIO,XL(*),XU(*),ACC,TNOW,PROFIL,PF(*)
      WRITE(IOUT,4)
      WRITE(IOUT,41) EVAL,TNOW,IFLAG  
      IF( O .EQ. 1 )THEN    
        IF(DABS(F).LE.1.0D+20)THEN
            WRITE(IOUT,42) F  
        ELSE
            WRITE(IOUT,82) F      
        ENDIF        
      ELSE   
        IF(DABS(F).LE.1.0D+20)THEN
            WRITE(IOUT,942) F  
        ELSE
            WRITE(IOUT,982) F      
        ENDIF 

        WRITE(IOUT,141) NINT(PF(1))

        DO I = 1,O
            IF(DABS(G(M+I)).LE.1.0D+14)THEN
              IF(G(M+O+I).LE.0.0D0)THEN
                WRITE(IOUT,70) I, G(M+I)
              ELSE
                WRITE(IOUT,70) I, -G(M+I) 
              ENDIF
            ELSE
              IF(G(M+O+I).LE.0.0D0)THEN
                WRITE(IOUT,71) I, G(M+I)
              ELSE
                WRITE(IOUT,71) I, -G(M+I) 
              ENDIF
            ENDIF
         ENDDO        
      ENDIF           
          IF(M.GT.0)THEN
            WRITE(IOUT,142)              
            IF(IFLAG.LT.100)THEN
                IF(VIO.LE.1.0D+12)THEN
                  WRITE(IOUT,183) VIO    
                ELSE
                  WRITE(IOUT,184) VIO           
                ENDIF
            ENDIF
          ENDIF
          IF(M.GT.0) WRITE(IOUT,144)              
          DO I = 1,ME
              IF(DABS(G(I)).LE.ACC)THEN
                  IF(DABS(G(I)).LE.1.0D+7)THEN
                      WRITE(IOUT,43) I,G(I)
                  ELSE
                      WRITE(IOUT,83) I,G(I)                  
                  ENDIF
              ELSE
                  IF(DABS(G(I)).LE.1.0D+7)THEN              
                      WRITE(IOUT,431) I,G(I)                  
                  ELSE
                      WRITE(IOUT,831) I,G(I)                  
                  ENDIF
              ENDIF
          ENDDO
          DO I = ME+1,M
              IF(G(I).GE.-ACC)THEN
                  IF(DABS(G(I)).LE.1.0D+7)THEN              
                      WRITE(IOUT,44) I,G(I)
                  ELSE
                      WRITE(IOUT,84) I,G(I)                  
                  ENDIF
              ELSE
                  IF(DABS(G(I)).LE.1.0D+7)THEN              
                      WRITE(IOUT,441) I,G(I)                  
                  ELSE
                      WRITE(IOUT,841) I,G(I)                  
                  ENDIF
              ENDIF
          ENDDO       
          WRITE(IOUT,145)                    
          DO I = 1,N
          PROFIL = -1.0D0
          IF( X(I).GT.XU(I)+1.0D-6 )THEN
              PROFIL = 91.0D0
              GOTO 88                  
          ENDIF    
          IF( X(I).LT.XL(I)-1.0D-6 )THEN
              PROFIL = 92.0D0
              GOTO 88                 
          ENDIF
          IF( XL(I).GT.XU(I) )THEN
              PROFIL = 93.0D0
              GOTO 88                  
          ENDIF                   
          IF( same(XL(I),XU(I)) .EQ. 1 )THEN
              PROFIL = 90.0D0
              GOTO 88                  
          ENDIF           
          IF( DABS(X(I)-XL(I)) .LE. (XU(I)-XL(I))/1000.0D0 )THEN
              PROFIL = 0.0D0
              GOTO 88
          ENDIF                
          IF( DABS(XU(I)-X(I)) .LE. (XU(I)-XL(I))/1000.0D0 )THEN
              PROFIL = 22.0D0
              GOTO 88                  
          ENDIF           
          DO J = 1,21      
              IF( X(I) .LE. XL(I) + DBLE(J) * (XU(I)-XL(I))/21.0D0)THEN
                  PROFIL = DBLE(J)
                  GOTO 88
              ENDIF
          ENDDO       
   88    CONTINUE
         IF(DABS(X(I)).LE.1.0D+14)THEN              
             IF(same(PROFIL, 0.0D0).EQ.1) WRITE(IOUT,400) I,X(I)
             IF(same(PROFIL, 1.0D0).EQ.1) WRITE(IOUT,401) I,X(I)
             IF(same(PROFIL, 2.0D0).EQ.1) WRITE(IOUT,402) I,X(I)
             IF(same(PROFIL, 3.0D0).EQ.1) WRITE(IOUT,403) I,X(I)
             IF(same(PROFIL, 4.0D0).EQ.1) WRITE(IOUT,404) I,X(I)
             IF(same(PROFIL, 5.0D0).EQ.1) WRITE(IOUT,405) I,X(I)
             IF(same(PROFIL, 6.0D0).EQ.1) WRITE(IOUT,406) I,X(I)
             IF(same(PROFIL, 7.0D0).EQ.1) WRITE(IOUT,407) I,X(I)
             IF(same(PROFIL, 8.0D0).EQ.1) WRITE(IOUT,408) I,X(I)
             IF(same(PROFIL, 9.0D0).EQ.1) WRITE(IOUT,409) I,X(I)
             IF(same(PROFIL,10.0D0).EQ.1) WRITE(IOUT,410) I,X(I)
             IF(same(PROFIL,11.0D0).EQ.1) WRITE(IOUT,411) I,X(I)
             IF(same(PROFIL,12.0D0).EQ.1) WRITE(IOUT,412) I,X(I)
             IF(same(PROFIL,13.0D0).EQ.1) WRITE(IOUT,413) I,X(I)
             IF(same(PROFIL,14.0D0).EQ.1) WRITE(IOUT,414) I,X(I)
             IF(same(PROFIL,15.0D0).EQ.1) WRITE(IOUT,415) I,X(I)
             IF(same(PROFIL,16.0D0).EQ.1) WRITE(IOUT,416) I,X(I)
             IF(same(PROFIL,17.0D0).EQ.1) WRITE(IOUT,417) I,X(I)
             IF(same(PROFIL,18.0D0).EQ.1) WRITE(IOUT,418) I,X(I)
             IF(same(PROFIL,19.0D0).EQ.1) WRITE(IOUT,419) I,X(I)
             IF(same(PROFIL,20.0D0).EQ.1) WRITE(IOUT,420) I,X(I)    
             IF(same(PROFIL,21.0D0).EQ.1) WRITE(IOUT,421) I,X(I)
             IF(same(PROFIL,22.0D0).EQ.1) WRITE(IOUT,422) I,X(I) 
             IF(same(PROFIL,90.0D0).EQ.1) WRITE(IOUT,490) I,X(I)
             IF(same(PROFIL,91.0D0).EQ.1) WRITE(IOUT,491) I,X(I)
             IF(same(PROFIL,92.0D0).EQ.1) WRITE(IOUT,492) I,X(I)
             IF(same(PROFIL,93.0D0).EQ.1) WRITE(IOUT,493) I,X(I)              
         ELSE      
             IF(same(PROFIL, 0.0D0).EQ.1) WRITE(IOUT,800) I,X(I)
             IF(same(PROFIL, 1.0D0).EQ.1) WRITE(IOUT,801) I,X(I)
             IF(same(PROFIL, 2.0D0).EQ.1) WRITE(IOUT,802) I,X(I)
             IF(same(PROFIL, 3.0D0).EQ.1) WRITE(IOUT,803) I,X(I)
             IF(same(PROFIL, 4.0D0).EQ.1) WRITE(IOUT,804) I,X(I)
             IF(same(PROFIL, 5.0D0).EQ.1) WRITE(IOUT,805) I,X(I)
             IF(same(PROFIL, 6.0D0).EQ.1) WRITE(IOUT,806) I,X(I)
             IF(same(PROFIL, 7.0D0).EQ.1) WRITE(IOUT,807) I,X(I)
             IF(same(PROFIL, 8.0D0).EQ.1) WRITE(IOUT,808) I,X(I)
             IF(same(PROFIL, 9.0D0).EQ.1) WRITE(IOUT,809) I,X(I)
             IF(same(PROFIL,10.0D0).EQ.1) WRITE(IOUT,810) I,X(I)
             IF(same(PROFIL,11.0D0).EQ.1) WRITE(IOUT,811) I,X(I)
             IF(same(PROFIL,12.0D0).EQ.1) WRITE(IOUT,812) I,X(I)
             IF(same(PROFIL,13.0D0).EQ.1) WRITE(IOUT,813) I,X(I)
             IF(same(PROFIL,14.0D0).EQ.1) WRITE(IOUT,814) I,X(I)
             IF(same(PROFIL,15.0D0).EQ.1) WRITE(IOUT,815) I,X(I)
             IF(same(PROFIL,16.0D0).EQ.1) WRITE(IOUT,816) I,X(I)
             IF(same(PROFIL,17.0D0).EQ.1) WRITE(IOUT,817) I,X(I)
             IF(same(PROFIL,18.0D0).EQ.1) WRITE(IOUT,818) I,X(I)
             IF(same(PROFIL,19.0D0).EQ.1) WRITE(IOUT,819) I,X(I)
             IF(same(PROFIL,20.0D0).EQ.1) WRITE(IOUT,820) I,X(I)    
             IF(same(PROFIL,21.0D0).EQ.1) WRITE(IOUT,821) I,X(I)
             IF(same(PROFIL,22.0D0).EQ.1) WRITE(IOUT,822) I,X(I) 
             IF(same(PROFIL,90.0D0).EQ.1) WRITE(IOUT,890) I,X(I)
             IF(same(PROFIL,91.0D0).EQ.1) WRITE(IOUT,891) I,X(I)
             IF(same(PROFIL,92.0D0).EQ.1) WRITE(IOUT,892) I,X(I)
             IF(same(PROFIL,93.0D0).EQ.1) WRITE(IOUT,893) I,X(I)             
          ENDIF         
      ENDDO   
      WRITE(IOUT,47)    
    4 FORMAT(' --------------------------------------------')
   41 FORMAT(' EVAL:',I9,',  TIME:',F9.2,',  IFLAG:',I4,/,
     &       ' --------------------------------------------')
   42 FORMAT(' F(X) =',F38.15)
   82 FORMAT(' F(X) =',E38.6) 
  942 FORMAT(' PROGRESS',F36.15)
  982 FORMAT(' PROGRESS',E36.6)      
  141 FORMAT(' --------------------------------------------',/,
     &       ' NUMBER OF PARETO POINTS',I21,/,
     &       ' --------------------------------------------')     
   70 FORMAT(' F(',I4,') =',F35.15)
   71 FORMAT(' F(',I4,') =',E35.8)       
  142 FORMAT(' --------------------------------------------')
  183 FORMAT(' VIOLATION OF G(X)',F27.12)  
  184 FORMAT(' VIOLATION OF G(X)',E27.6)  
  144 FORMAT(' --------------------------------------------')   
   43 FORMAT(' G(',I4,') =',F16.8,'  (EQUALITY CONSTR)')
   44 FORMAT(' G(',I4,') =',F16.8,'  (IN-EQUAL CONSTR)')
  431 FORMAT(' G(',I4,') =',F16.8,'  (EQUALITY CONSTR)',
     &       '  <---  INFEASIBLE  ( G NOT = 0 )')
  441 FORMAT(' G(',I4,') =',F16.8,'  (IN-EQUAL CONSTR)',
     &       '  <---  INFEASIBLE  ( G < 0 )')
   83 FORMAT(' G(',I4,') =',E16.3,'  (EQUALITY CONSTR)')
   84 FORMAT(' G(',I4,') =',E16.3,'  (IN-EQUAL CONSTR)')
  831 FORMAT(' G(',I4,') =',E16.3,'  (EQUALITY CONSTR)',
     &       '  <---  INFEASIBLE  ( G NOT = 0 )')
  841 FORMAT(' G(',I4,') =',E16.3,'  (IN-EQUAL CONSTR)',
     &       '  <---  INFEASIBLE  ( G < 0 )')     
  145 FORMAT(' --------------------------------------------',
     &                              '            BOUNDS-PROFIL    ')
  400 FORMAT(' X(',I5,') = ',F33.15,'    !   XL___________________')
  401 FORMAT(' X(',I5,') = ',F33.15,'    !   x____________________')
  402 FORMAT(' X(',I5,') = ',F33.15,'    !   _x___________________')
  403 FORMAT(' X(',I5,') = ',F33.15,'    !   __x__________________')
  404 FORMAT(' X(',I5,') = ',F33.15,'    !   ___x_________________')
  405 FORMAT(' X(',I5,') = ',F33.15,'    !   ____x________________')
  406 FORMAT(' X(',I5,') = ',F33.15,'    !   _____x_______________')
  407 FORMAT(' X(',I5,') = ',F33.15,'    !   ______x______________')
  408 FORMAT(' X(',I5,') = ',F33.15,'    !   _______x_____________')
  409 FORMAT(' X(',I5,') = ',F33.15,'    !   ________x____________')
  410 FORMAT(' X(',I5,') = ',F33.15,'    !   _________x___________')
  411 FORMAT(' X(',I5,') = ',F33.15,'    !   __________x__________')
  412 FORMAT(' X(',I5,') = ',F33.15,'    !   ___________x_________')
  413 FORMAT(' X(',I5,') = ',F33.15,'    !   ____________x________')
  414 FORMAT(' X(',I5,') = ',F33.15,'    !   _____________x_______')
  415 FORMAT(' X(',I5,') = ',F33.15,'    !   ______________x______')
  416 FORMAT(' X(',I5,') = ',F33.15,'    !   _______________x_____')
  417 FORMAT(' X(',I5,') = ',F33.15,'    !   ________________x____')
  418 FORMAT(' X(',I5,') = ',F33.15,'    !   _________________x___')
  419 FORMAT(' X(',I5,') = ',F33.15,'    !   __________________x__')
  420 FORMAT(' X(',I5,') = ',F33.15,'    !   ___________________x_')
  421 FORMAT(' X(',I5,') = ',F33.15,'    !   ____________________x')  
  422 FORMAT(' X(',I5,') = ',F33.15,'    !   ___________________XU')  
  490 FORMAT(' X(',I5,') = ',F33.15,'    !   WARNING: XL = XU     ')  
  491 FORMAT(' X(',I5,') = ',F33.15,'  <---  *** ERROR *** (X > XU) ')  
  492 FORMAT(' X(',I5,') = ',F33.15,'  <---  *** ERROR *** (X < XL) ')  
  493 FORMAT(' X(',I5,') = ',F33.15,'  <---  *** ERROR *** (XL > XU)')
  800 FORMAT(' X(',I5,') = ',E33.1,'    !   XL___________________')
  801 FORMAT(' X(',I5,') = ',E33.1,'    !   x____________________')
  802 FORMAT(' X(',I5,') = ',E33.1,'    !   _x___________________')
  803 FORMAT(' X(',I5,') = ',E33.1,'    !   __x__________________')
  804 FORMAT(' X(',I5,') = ',E33.1,'    !   ___x_________________')
  805 FORMAT(' X(',I5,') = ',E33.1,'    !   ____x________________')
  806 FORMAT(' X(',I5,') = ',E33.1,'    !   _____x_______________')
  807 FORMAT(' X(',I5,') = ',E33.1,'    !   ______x______________')
  808 FORMAT(' X(',I5,') = ',E33.1,'    !   _______x_____________')
  809 FORMAT(' X(',I5,') = ',E33.1,'    !   ________x____________')
  810 FORMAT(' X(',I5,') = ',E33.1,'    !   _________x___________')
  811 FORMAT(' X(',I5,') = ',E33.1,'    !   __________x__________')
  812 FORMAT(' X(',I5,') = ',E33.1,'    !   ___________x_________')
  813 FORMAT(' X(',I5,') = ',E33.1,'    !   ____________x________')
  814 FORMAT(' X(',I5,') = ',E33.1,'    !   _____________x_______')
  815 FORMAT(' X(',I5,') = ',E33.1,'    !   ______________x______')
  816 FORMAT(' X(',I5,') = ',E33.1,'    !   _______________x_____')
  817 FORMAT(' X(',I5,') = ',E33.1,'    !   ________________x____')
  818 FORMAT(' X(',I5,') = ',E33.1,'    !   _________________x___')
  819 FORMAT(' X(',I5,') = ',E33.1,'    !   __________________x__')
  820 FORMAT(' X(',I5,') = ',E33.1,'    !   ___________________x_')
  821 FORMAT(' X(',I5,') = ',E33.1,'    !   ____________________x')  
  822 FORMAT(' X(',I5,') = ',E33.1,'    !   ___________________XU')  
  890 FORMAT(' X(',I5,') = ',E33.1,'    !   WARNING: XL = XU     ')  
  891 FORMAT(' X(',I5,') = ',E33.1,'  <---  *** ERROR *** (X > XU) ')  
  892 FORMAT(' X(',I5,') = ',E33.1,'  <---  *** ERROR *** (X < XL) ')  
  893 FORMAT(' X(',I5,') = ',E33.1,'  <---  *** ERROR *** (XL > XU)')   
   47 FORMAT(/,' ')    
      END
      SUBROUTINE PRINT_FINAL(IFLAG,TNOW,TMAX,EVAL,MAXEVAL,O,N,M,ME,
     &                       X,G,F,XL,XU,RW,ACC,WRES,PF,PARAM,IOUT) 
      IMPLICIT NONE
      INTEGER IFLAG,EVAL,MAXEVAL,O,N,M,ME,WRES,IOUT
      DOUBLE PRECISION TNOW,TMAX,X(*),G(*),F,XL(*),XU(*),RW(*),ACC
      DOUBLE PRECISION PARAM(*),PF(*)
      IF(IFLAG.EQ.1.OR.IFLAG.EQ.2)THEN 
        IF(TNOW.GE.TMAX)    WRITE(IOUT,411)  
        IF(EVAL.GE.MAXEVAL) WRITE(IOUT,412)
      ENDIF
      IF(IFLAG.EQ.3.OR.IFLAG.EQ.4) WRITE(IOUT,413) NINT(PARAM(4))
      IF(IFLAG.EQ.5.OR.IFLAG.EQ.6) WRITE(IOUT,414) NINT(PARAM(5))
      IF(IFLAG.EQ.7) WRITE(IOUT,415)      
 411  FORMAT(/,' OPTIMIZATION FINISHED  --->  MAXTIME REACHED')
 412  FORMAT(/,' OPTIMIZATION FINISHED  --->  MAXEVAL REACHED')
 413  FORMAT(/,' OPTIMIZATION FINISHED  --->  ALGOSTOP (=',I3,')')
 414  FORMAT(/,' OPTIMIZATION FINISHED  --->  EVALSTOP (=',I9,')')
 415  FORMAT(/,' OPTIMIZATION FINISHED  --->  FSTOP REACHED')
      WRITE(IOUT,42)    
   42 FORMAT(/,/,'         BEST SOLUTION FOUND BY MIDACO       ')       
      CALL PRINT_SOLUTION( O, N, M, ME, X, G, F, RW(WRES), PF,
     &                     XL, XU, ACC, EVAL, TNOW, IFLAG, IOUT)
      END
      SUBROUTINE WARNINGS_AND_ERRORS( IFLAG , IOUT )     
      IMPLICIT NONE
      INTEGER IFLAG,IOUT
      IF(IFLAG.LT.100)THEN
      WRITE(IOUT,1) IFLAG
    1 FORMAT(/,' *** WARNING ***   ( IFLAG =',I6,' )',/)              
      ELSE
      WRITE(IOUT,2) IFLAG
    2 FORMAT(/,/,/,' *** MIDACO INPUT ERROR ***   ( IFLAG =',I6,' )')   
      ENDIF             
      END
      SUBROUTINE PRINT_PARETOFRONT(O,M,N,PF,BESTX,BESTG,PFMAX,IOUT)
      IMPLICIT NONE
      INTEGER O,M,N,  IOUT,K,I,PSIZE,PFMAX
      DOUBLE PRECISION PF(*),BESTX(*),BESTG(*),DUMMY
      PSIZE = INT(PF(1))     
      OPEN(IOUT,FILE='MIDACO_PARETOFRONT.tmp',STATUS='UNKNOWN')  
      WRITE(IOUT,1)
    1 FORMAT(
     & '#########################################################',/,
     & '### This file contains the pareto front approximation ###',/,
     & '#########################################################',/,
     & '### Solution format:     F(1:O)    G(1:M)    X(1:N)   ###',/,
     & '#########################################################')
      write(IOUT,2) O,m,n,int(pf(1))
    2 format('#',/,
     &       '#        O         M         N     PSIZE',/,
     &       '#',/,
     &       4I10)  
      write(IOUT,7)
    7 format('#',/,
     &       '#        MIDACO solution',/,
     &       '#')
        do i=1,O
          
          dummy = bestg(m+i)
          if( bestg(m+O+i) .GT. 0.0D0 ) dummy = -dummy

          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy           
        enddo
        do i=1,M
          dummy =  bestg(i)
          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy 
        enddo
        do i=1,N
          dummy =  bestx(i)
          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy 
        enddo       
      write(IOUT,3)
    3 format(/,'#',/,
     &         '#        All non-dominated solutions found by MIDACO',/,
     &         '#')
      do k=1,PSIZE 
        do i=1,O
          dummy = pf( 2 + O*(k-1)+i-1 )
          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy           
        enddo
        do i=1,M
          dummy =  pf( 2 + O*PFMAX + M*(k-1)+i-1 ) 
          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy 
        enddo
        do i=1,N
          dummy =  pf( 2 + O*PFMAX + M*PFMAX + N*(k-1)+i-1 ) 
          if( DABS(dummy) .LT. 1.0D9) write(IOUT, 4) dummy
          if( DABS(dummy) .GE. 1.0D9) write(IOUT,44) dummy 
        enddo        
        write(IOUT,5)
      enddo
    4 format(F19.7,$) 
   44 format(E16.5,$) 
    5 format(' ')    
      CALL CLOSE_OUTPUT(IOUT)
      call rename('MIDACO_PARETOFRONT.tmp','MIDACO_PARETOFRONT.TXT')  
      END
      SUBROUTINE SAVE_HISTORY(SAVE2FILE,P,O,N,M,F,G,X,IOUT,IFLAG,ISTOP)     
      INTEGER SAVE2FILE,P,O,N,M,IOUT ,I,K ,TIX,IFLAG,ISTOP
      DOUBLE PRECISION F(*),G(*),X(*),PREVIOUS_X(10000),DUMMY     
      DATA PREVIOUS_X /10000 * 0.0d0/
      DATA TIX /0/
      IF(IFLAG.EQ.0) TIX = 0
      TIX = TIX + 1
      IF(TIX.GT.SAVE2FILE) TIX = 1
      IF(TIX.EQ.1)THEN
C       Open and Delete file, otherwise delete-only might crash        
        OPEN(IOUT,FILE='MIDACO_HISTORY.TXT',STATUS='UNKNOWN')
        CLOSE(IOUT,STATUS='DELETE')
      ENDIF      
      IF(TIX.EQ.1)THEN
      OPEN(IOUT,FILE='MIDACO_HISTORY.TXT'  ,STATUS='UNKNOWN')
      WRITE(IOUT,1)
    1 FORMAT(
     & '###################################################',/,
     & '### This file contains the history of solutions ###',/,
     & '###################################################',/,
     & '### SOLUTION FORMAT:  F(1:O)  G(1:M)  X(1:N)    ###',/,
     & '###################################################')
      WRITE(IOUT,2) O,M,N
    2 FORMAT('#',/,
     &       '#        O         M         N  ',/,
     &       '#',/,
     &       3I10)  
      WRITE(IOUT,3)
    3 FORMAT(/,'#',/,
     &         '#        SOLUTION HISTORY (in chronological order)',/,
     &         '#')   
      ELSE
      IF(TIX.EQ.2)THEN
             DO I=1,O
               DUMMY = F(I)
               IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
               IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY           
             ENDDO
             DO I=1,M
               DUMMY =  G(I)
               IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
               IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY 
             ENDDO       
             DO I=1,N
               DUMMY =  PREVIOUS_X(I)
               IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
               IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY 
             ENDDO 
             WRITE(IOUT,5)  
      ELSE
      DO K=1,P
        DO I=1,O
          IF(ISTOP.LE.0) DUMMY = F((K-1)*O+I)
          IF(ISTOP.GE.1) DUMMY = F(I)
          IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
          IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY           
        ENDDO
        DO I=1,M
          IF(ISTOP.LE.0) DUMMY = G((K-1)*M+I)
          IF(ISTOP.GE.1) DUMMY = G(I)
          IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
          IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY 
        ENDDO       
        DO I=1,N
          IF(ISTOP.LE.0) DUMMY = PREVIOUS_X((K-1)*N+I)
          IF(ISTOP.GE.1) DUMMY = X(I)
          IF( DABS(DUMMY) .LT. 1.0D9) WRITE(IOUT, 4) DUMMY
          IF( DABS(DUMMY) .GE. 1.0D9) WRITE(IOUT,44) DUMMY 
        ENDDO
        WRITE(IOUT,5)
      ENDDO
    4 FORMAT(F19.7,$) 
   44 FORMAT(E16.5,$) 
    5 FORMAT(' ') 
      ENDIF
      ENDIF
      DO K=1,P
       DO I=1,N
        PREVIOUS_X((K-1)*N+I) = X((K-1)*N+I)
       ENDDO   
      ENDDO  
      CALL FLUSH_OUTPUT( IOUT )      
      IF(ISTOP.GE.1) CALL CLOSE_OUTPUT( IOUT ) 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      integer function same( a, b )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision a, b
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      same = 0
      if( a.lt.b ) return
      if( a.gt.b ) return      
      same = 1      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine calculates the cpu-time.
C     In case of compiler problems, you may want to change this routine.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GET_TIME( SECONDS, PARALLEL, PTYPE )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION SECONDS
      INTEGER PARALLEL, PTYPE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL CPU_TIME( SECONDS )         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     For parallelization via openMP, the time command is adjusted
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF( PTYPE .EQ. 1 ) SECONDS = SECONDS / DBLE(PARALLEL)   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CLOSE_OUTPUT( IOUT )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER IOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CLOSE( IOUT )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine forces to flush the MIDACO output to text files.
c     Some Fortran compiler (e.g. g77 and f77) do not accept the 'flush'
c     command. In case you have any problems with 'flush', you can
c     savely REMOVE or COMMENT-OUT the below 'flush' command.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLUSH_OUTPUT( IOUT )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER IOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FLUSH( IOUT )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END
c     END OF FILE                

  