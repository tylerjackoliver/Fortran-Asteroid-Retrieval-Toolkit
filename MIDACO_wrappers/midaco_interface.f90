module midaco_interface

    implicit none

    interface

        INTEGER (C_INT) function midaco_c(p, o, n, ni, m, me, x, f, g, xl, xu, iflag, &
                            istop, param, rw, lrw, iw, liw, pf, lpf, save2file, maxeval, &
                            maxtime, printeval) bind(C, name="midaco_wrap")

            use iso_c_binding

            implicit none

            INTEGER(kind=C_LONG)        :: p            ! Parallelisation factor
            INTEGER(kind=C_LONG)        :: o            ! Number of objective function
            INTEGER(kind=C_LONG)        :: n            ! Number of optimisation variable
            INTEGER(kind=C_LONG)        :: ni           ! Number of integer optimisation variables
            INTEGER(kind=C_LONG)        :: m            ! Number of constraints in total
            INTEGER(kind=C_LONG)        :: me           ! Number of equality constraints

            TYPE(C_PTR), VALUE          :: x            ! Array containing the n-th iterate
            TYPE(C_PTR), VALUE          :: f            ! Array containing the objective function values
            TYPE(C_PTR), VALUE          :: g            ! Array containing constraint values G
            TYPE(C_PTR), VALUE          :: xl           ! Array containing lower bounds
            TYPE(C_PTR), VALUE          :: xu           ! Array containing upper bounds
            
            INTEGER(kind=C_LONG)        :: iflag        ! Information flag used by MIDACO
            INTEGER(kind=C_LONG)        :: istop        ! Communication flag used by MIDACO
            
            TYPE(C_PTR), VALUE          :: param        ! Array containing 13 parameters that can be selected by the user
            TYPE(C_PTR), VALUE          :: rw           ! Real workarray
            
            INTEGER(kind=C_LONG)        :: lrw          ! Length of the real workarray
            
            TYPE(C_PTR), VALUE          :: iw           ! Integer workarray
            
            INTEGER(kind=C_LONG)        :: liw          ! Length of integer workarray
        
            TYPE(C_PTR), VALUE          :: pf           ! Pareto front approximation

            INTEGER(kind=C_LONG)        :: lpf          ! Length of the pareto front
            INTEGER(kind=C_LONG)        :: save2file    ! 
            INTEGER(kind=C_LONG)        :: maxeval  
            INTEGER(kind=C_LONG)        :: maxtime
            INTEGER(kind=C_LONG)        :: printeval

        end function midaco_c

        INTEGER (C_INT) function reset_midaco_run() bind(C, name="reset_run")

            use iso_c_binding
            implicit none

        end function reset_midaco_run

    end interface

    contains

        subroutine midaco(p, o, n, ni, m, me, x, f, g, xl, xu, iflag, &
                    istop, param, rw, lrw, iw, liw, pf, lpf, save2file, &
                    maxeval, maxtime, printeval)

            use iso_c_binding

            implicit none

            INTEGER, INTENT(IN)                         :: p                                ! Parallelisation factor
            INTEGER, INTENT(IN)                         :: o                                ! Number of objective function
            INTEGER, INTENT(IN)                         :: n                                ! Number of optimisation variable
            INTEGER, INTENT(IN)                         :: ni                               ! Number of integer optimisation variables
            INTEGER, INTENT(IN)                         :: m                                ! Number of constraints in total
            INTEGER, INTENT(IN)                         :: me                               ! Number of equality constraints

            DOUBLE PRECISION, INTENT(INOUT)             :: x(:)                             ! Array containing the n-th iterate
            DOUBLE PRECISION, INTENT(INOUT)             :: f(:)                             ! Array containing the objective function values
            DOUBLE PRECISION, INTENT(INOUT)             :: g(:)                             ! Array containing constraint values G
            DOUBLE PRECISION, INTENT(INOUT)             :: xl(:)                            ! Array containing lower bounds
            DOUBLE PRECISION, INTENT(INOUT)             :: xu(:)                            ! Array containing upper bounds
            
            INTEGER, INTENT(INOUT)                      :: iflag                            ! Information flag used by MIDACO
            INTEGER, INTENT(INOUT)                      :: istop                            ! Communication flag used by MIDACO
            
            DOUBLE PRECISION, INTENT(INOUT)             :: param(:)                         ! Array containing 13 parameters that can be selected by the user
            DOUBLE PRECISION, INTENT(INOUT)             :: rw(:)                            ! Real workarray
            
            INTEGER, INTENT(IN)                         :: lrw                              ! Length of the real workarray
            
            INTEGER, INTENT(INOUT)                      :: iw(:)                            ! Integer workarray
            
            INTEGER, INTENT(IN)                         :: liw                              ! Length of integer workarray
        
            DOUBLE PRECISION, INTENT(INOUT)             :: pf(:)                            ! Pareto front approximation

            INTEGER, INTENT(IN)                         :: lpf                              ! Length of the pareto front

            INTEGER, INTENT(IN)                         :: save2file                        ! Determines result behaviour
            INTEGER, INTENT(IN)                         :: maxeval                          ! Maximum function evaluations
            INTEGER, INTENT(IN)                         :: maxtime                          ! Maximum evaluation time
            INTEGER, INTENT(IN)                         :: printeval

            INTEGER(KIND=C_INT)                         :: RESULT

            INTEGER(kind=C_LONG)                        :: pISO                             ! Parallelisation factor
            INTEGER(kind=C_LONG)                        :: oISO                             ! Number of objective function
            INTEGER(kind=C_LONG)                        :: nISO                             ! Number of optimisation variable
            INTEGER(kind=C_LONG)                        :: niISO                            ! Number of integer optimisation variables
            INTEGER(kind=C_LONG)                        :: mISO                             ! Number of constraints in total
            INTEGER(kind=C_LONG)                        :: meISO                            ! Number of equality constraints

            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: xISO(:)                          ! Array containing the n-th iterate
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: fISO(:)                          ! Array containing the objective function values
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: gISO(:)                          ! Array containing constraint values G
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: xlISO(:)                         ! Array containing lower bounds
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: xuISO(:)                         ! Array containing upper bounds
            
            INTEGER(kind=C_LONG)                        :: iflagISO                         ! Information flag used by MIDACO
            INTEGER(kind=C_LONG)                        :: istopISO                         ! Communication flag used by MIDACO
            
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: paramISO(:)                      ! Array containing 13 parameters that can be selected by the user
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: rwISO(:)                         ! Real workarray
            
            INTEGER(kind=C_LONG)                        :: lrwISO                           ! Length of the real workarray
            
            INTEGER(kind=C_LONG), ALLOCATABLE, TARGET   :: iwISO(:)                         ! Integer workarray
            
            INTEGER(kind=C_LONG)                        :: liwISO                           ! Length of integer workarray
        
            REAL(kind=C_DOUBLE), ALLOCATABLE, TARGET    :: pfISO(:)                         ! Pareto front approximation

            INTEGER(kind=C_LONG)                        :: lpfISO                           ! Length of the pareto front

            INTEGER(kind=C_LONG)                        :: save2fileISO
            INTEGER(kind=C_LONG)                        :: maxevalISO
            INTEGER(kind=C_LONG)                        :: maxtimeISO
            INTEGER(kind=C_LONG)                        :: printevalISO

            ! Allocate ALLOCATABLE, TARGETsed on input array            
            
            ALLOCATE(xISO(size(x)),         &
                     fISO(size(f)),         &
                     gISO(size(g)),         &
                     xlISO(size(xl)),       & 
                     xuISO(size(xu)),       &
                     paramISO(size(param)), &
                     rwISO(size(rw)),       &
                     iwISO(size(iw)),       &
                     pfISO(size(pf)))

            ! Convert reals

            xISO = x; fISO = f; gISO = g; xlISO = xl; xuISO = xu; paramISO = param; rwISO = rw; iwISO = iw; pfISO = pf;

            ! Convert integer kinds

            pISO            = INT(p, C_LONG)
            oISO            = INT(o, C_LONG)
            nISO            = INT(n, C_LONG)
            niISO           = INT(ni, C_LONG)
            mISO            = INT(m, C_LONG)
            meISO           = INT(me, C_LONG)
            iflagISO        = INT(iflag, C_LONG)
            lpfISO          = INT(lpf, C_LONG)
            istopISO        = INT(istop, C_LONG)
            lrwISO          = INT(lrw, C_LONG)
            liwISO          = INT(liw, C_LONG)
            save2fileISO    = INT(save2file, C_LONG)
            maxevalISO      = INT(maxeval, C_LONG)
            maxtimeISO      = INT(maxtime, C_LONG)
            printevalISO    = INT(printeval, C_LONG)

            ! Call the MIDACO routines
            ! C_LOC gets the C pointer to the address

            RESULT =  MIDACO_c(pISO, oISO, nISO, niISO, mISO, meISO, C_LOC(xISO), C_LOC(fISO), C_LOC(gISO), &
                               C_LOC(xlISO), C_LOC(xuISO), iflagISO, istopISO, C_LOC(paramISO), C_LOC(rwISO), &
                               lrwISO, C_LOC(iwISO), liwISO, C_LOC(pfISO), lpfISO, save2fileISO, maxevalISO, &
                               maxtimeISO, printevalISO)

            ! Now convert back

            x = xISO; f = fISO; g = gISO; xl = xlISO; xu = xuISO; param = paramISO; rw = rwISO; iw = iwISO; pf = pfISO;
            iflag = INT(iflagISO); istop = INT(istopISO); 

            ! Deallocate ALLOCATABLE, TARGET           

            DEALLOCATE(xISO,        &
                       fISO,        &
                       gISO,        &
                       xlISO,       & 
                       xuISO,       &
                       paramISO,    &
                       rwISO,       &
                       iwISO,       &
                       pfISO)   

        end subroutine midaco
        

end module midaco_interface