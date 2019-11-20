! ============================================================
!
!
!
!     ASTEROID PRUNING PRE-FILTER: FORTRAN08 IMPLEMENTATION
!
!
!
! =============================================================
!
!   INTRODUCTION
!   ------------
!   The below program is essentially four nested do-loops that apply the 
!   pre-filtering methodology to every candidate in the minor bodies
!   database and target OE files.
!
!   LIST OF DEPENDENCIES
!   --------------------
!   OpenMP
!   gFortran with --free-fform
!   Minor bodies database as .csv
!   Target OE files as .csv
!
!   METADATA
!   --------
!   Date: 31/03/18
!   Author: Jack Tyler
!
!   CHANGELOG
!   ---------
!   1/1/18: Bug fixes
!   2/2/2019: Minor refactoring
!


program main

    use omp_lib
    use constants                                                                       ! Standardises constants for use in all Fortran codes
    use iso_fortran_env, only : error_unit 

    implicit none

    integer                    :: nthreads = 2                                          ! Number of threads to use
    integer                    :: unitno, candidatecount, itercount                     ! Counting variables
    integer                    :: i, j, k                                               ! Loop variables 
    integer                    :: bodyID                                                ! SPK ID of the body being investigated
    integer                    :: desireddate                                           ! Week ID to compute the solution for
                                                                                        ! For use on Lyceum, the input was parameterised by date. At higher core counts,
                                                                                        ! It becomes feasible to do all dates at once.

    real                       :: a_b                                                   ! SMA of body
    real                       :: e_b                                                   ! Eccentricity of body 
    real                       :: inb                                                   ! Inclination of body (rad!!)
    real                       :: rp_t                                                  ! Periapse radius of target
    real                       :: e_t
    real                       :: in_t                                                  ! Inclination of target (deg!!)
    real                       :: start                                                 ! OMP WTIME timing
    real                       :: finish                                                ! ""
    real                       :: M_t                                                   ! Mean anomaly of the target
    real                       :: mindV                                                 ! minimum delta V
    real                       :: t1, t2, t3, t4, t5, date, o, Om, M, p                 ! Temporary variables


    integer, dimension(17719)  :: database                                              ! SPK ID Database
    ! integer, dimension(32) :: database

    real, dimension(1353, 3)   :: bodyData                                              ! rp_b, e_b, inb,  of the target body (needed only for Hohmann)
    real, dimension(360000,4)  :: targetData                                            ! rp_b, e_t, in_t, M_t of the target manifold section

    character(len=260)         :: bodystring, filepath                                  ! String of the body SPK ID (file IO), filepath for directory prefix
    character(len=8)           :: fmt = '(I7.7)'                                        ! Formatting for result I/O

    ! Variable initialisations

    candidatecount = 0                                                                  ! Number of candidates found                                                          
    itercount = 0                                                                       ! Number of weeks processed
    filepath = 'corrected_oe/'                                                          ! Directory prefix: alter as req'd

    ! Print welcome message

    write(*,*) "============================================================"
    write(*,*) " "
    write(*,*) "ASTEROID PRUNING PRE-FILTER: FORTRAN08 IMPLEMENTATION"
    write(*,*) " "
    write(*,*) "LYCEUM VERSION: LIMITED PRINTING CAPABILITIES AVAILABLE"
    write(*,*) "DESIGNED FOR USE WITH 16 THREADS"
    write(*,*) " "
    write(*,*) "============================================================"

    ! Initialise OMP and delete current results.txt

!    call OMP_SET_NUM_THREADS(nthreads)    
!   call file_check()
    ! call getarg(1, desireddate)                                                       ! If parameterising by input date: get date of interest

    ! Load data sets

    write(*,*) "Loading data sets..."

    open(unit=3, file="L2_endConds.csv")                                                 ! Target OEs

    do i = 1,360000

        read(3,*) targetData(i,1), targetData(i,2), targetData(i,3), t1, t2, targetData(i,4), t3, t4

    end do

    close(3)

    open(unit=20, file="Database_SPKs.csv")                                             ! Database file
    
    read(20,*)                                                                          ! Read once to avoid 'spk_id' header

    do i = 1,size(database)
    
        read(20,*) database(i)
    
    end do

    close(20)

    ! Write progress message and begin timing

    write(*,*) "Loading complete. Beginning compute sequence for body #", desireddate
!    start = omp_get_wtime()

    !$OMP PARALLEL DO PRIVATE(I,J,K, bodyid, bodystring, unitno, mindV, bodyData, a_b, e_b, inb, rp_t, e_t, in_t)
    do i = 1, size(database) ! All bodies

        ! Get the body ID from the database and write to a string

        bodyID = database(i)
        write(bodystring, fmt) bodyID

        ! Unique unit number for each thread - need to do this every loop?

 !       unitno = OMP_GET_THREAD_NUM() + 30
        unitno = 30
        ! Open file, read the input data and close again

        open(unit=unitno, file=trim(filepath)//trim(bodystring)//'_OE_corr.csv')
    
        do j = 1,1353

            read(unitno,*) date, bodyData(j,1), bodyData(j,2), bodyData(j,3), o, O, M, p

        end do

        close(unitno)

        ! Loop at the desired data through every body, and every target, to see if we have any possible transfers

        ! exitloop: do j = 272,1353
        exitloop:  do j = 221,221

            a_b = bodyData(j,1)*au                                                                                  ! a is non-dim
            e_b = bodyData(j,2)
            inb = bodyData(j,3)                                                                                     ! Note: degrees - is converted in cheapest subroutine

            do k = 1, 360000
                
                rp_t = targetData(k,1)*au                                                                           ! rp_t is non-dim
                e_t = targetData(k,2)
                in_t = targetData(k,3)
                M_t = targetData(k,4)

                call cheapest(a_b, e_b, inb, rp_t, e_t, in_t, M_t, mindV)                                           ! Compute if we have any possible transfers

                if (mindV < .7) then 

                    call file_add(bodyID, unitno, mindV, k)                                                         ! Open results file, add the SPK ID in, and close again
                    candidatecount = candidatecount + 1
		            exit exitloop                                                                                   ! If we have a candidate, then we don't care about any other transfer opportunities; stop.

                end if

            end do

        end do exitloop

        itercount = itercount + 1

        write(*,*) "Completed processing object number ", itercount
        write(*,*) "Candidates found: ", candidatecount

    end do
    !$OMP END PARALLEL DO

!    finish = omp_get_wtime()

    write(*,*) "==================== Scan completed ========================"
    write(*,*) "Objects scanned: ", itercount
    write(*,*) "Candidates found: ", candidatecount
!    write(*,*) "Compute time: ", finish-start, "seconds."

end program main


!
! Routine to add the SPK ID of a body to a file if it is deemed to be a candidate
! WARNING: May have issues if multiple threads attempt to I/O at once.
!

subroutine file_add(bodyID, unitnum, mindv, k)

    ! Subroutine intents

    character(len=13) :: filepath ! Define address of file to be read
    integer, intent(in) :: bodyID
    integer, intent(in) :: unitnum
    real, intent(in) :: mindv
    integer, intent(in) :: k
    integer :: outval1, outval2

    character(len=11) :: outstring

    logical :: exist

    ! call where_from(k, outstring)

    filepath = "results.txt"

    inquire(file=filepath, exist=exist)

    if (exist) then

      open(unitnum, file=filepath, status="old", position="append", action="write")

    else

      open(unitnum, file=filepath, status="new", action="write")

    end if

    ! Get upper and lower bound of where the orbit family is (hardcoded bounds from generation)

    call data_index(k, outval1, outval2)

    ! Write the body, the mindv, and the family (bounds)

    write(unitnum, *) bodyID, mindv, outval1, outval2
    close(unitnum)

end subroutine file_add

!
! Subroutine to check if the results file already exists in the working folder
!

subroutine file_check()

    character(len=13) :: filepath
    logical :: exist

    filepath = 'results.txt'

    inquire(file=filepath, exist=exist)

    if (exist) then

        open(20, file=filepath)
        close(20, status='delete')

    end if

end subroutine file_check


subroutine cheapest(a_b, e_b, inb, rp_t, e_t, in_t, M_t, cheap)

    real, intent(in)    :: a_b                  ! SMA of the candidate (km)
    real, intent(in)    :: e_b                  ! Eccentricity of the candidate
    real, intent(in)    :: inb                  ! Inclination of the candiate (rad)
    real, intent(in)    :: rp_t                 ! Radius of periapsis of the target (km)
    real, intent(in)    :: e_t                  ! Eccentricity of the target
    real, intent(in)    :: in_t                 ! Inclination of the target (rad)
    real, intent(in)    :: M_t                  ! Mean anomaly of the target (rad)

    real, intent(out)   :: cheap                ! Cheapest Hohmann transfer (km/s)

    real                :: mu_s                 ! Solar graviational parameter (km ^ 2/s ^ 3)
    real                :: ra_b                 ! Radius of apoapsis of candidate (km)
    real                :: rp_b                 ! Radius of periapsis of candidate (km)
    real                :: a_int                ! SMA of intermediate transfer arc (km)
    real                :: a_f                  ! Equivalent SMA of target (doesn't really exist)
    real                :: rstar                ! Quotient for inclination change calculation
    real                :: dVi                  ! Velocity required for inclination change (km/s)
    real                :: dV1                  ! Velocity required for first Hohmann maneouvre (km/s)
    real                :: dV2                  ! Velocity required for second Hohmann maneouvre (km/s)
    real                :: v1                   ! Total velocity required for first maneouvre (km/s)
    real                :: v2                   ! Total velocity required for second maneouvre (km/s)
    real                :: rt                   ! Perifocal radius of the target (km)
    real                :: theta_t              ! True anomaly of the target (rad)

    double precision    :: pi = 4.*datan(1.d0)

    mu_s = 1.327124400189e11

    ! Convert inb to radian; print warning message if (poorly-constructed) test for degrees fails.

!    if (abs(inb) < 2*pi) then

!        write(error_unit,*) "Warning: inb should be passed in as degrees; inb might be in radians here."  ! Writes to stderr on gfortran

!   end if

    in_b = inb*pi/180.0

    !
    ! We need to convert the mean anomaly of the target into the true anomaly
    ! to get the perifocal distance
    !

    call mean_to_true(M_t, e, theta_t)

    ! Now get the perifocal distance of the target

    r_t = rp_t * (1 + e) / (1 + e_t*cos(theta_t))

    ! First, Compute apoapsis and periapsis distance for the body
    
    ra_b = a_b * (1 + e_b)
    rp_b = a_b * (1 - e_b)

    !
    ! There are four possible cases of transfers:
    !       1. Ap of body -> target, inclination change at body
    !       2. Ap of body -> target, inclination change at target
    !       3. Peri of body -> target, inclination change at body
    !       4. Peri of body -> target, inclination change at target
    !
    ! We can, however, automatically eliminate two of these cases by simply choosing the most efficient point
    ! to do the inclination change.
    !

    !
    ! Case 1/2
    !

    ! First, compute the intermediate and final semimajor axes

    a_0 = a_b
    a_int = (ra_b + r_t)
    a_f = rp_t / (1 - e_t)

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/ra_b - 1/a_int)) - sqrt(mu_s * (2/ra_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/r_t - 1/a_f)) - sqrt(mu_s * (2/r_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (ra_b > r_t) then  ! If ra_b is larger, then it's better to do the inclination change post-apside change

        rstar = r_t / ra_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b))

        v1 = sqrt(dV1 ** 2 + dVi **2) + abs(dV2)

    else

        rstar = ra_b/r_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b))
        v1 = abs(dV1) + sqrt(dV2 ** 2 + dVi **2)

    end if

    !
    ! Case 3/4
    !

    ! Update the intermediate SMA

    a_int = (rp_b + r_t)

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/rp_b - 1/a_int)) - sqrt(mu_s * (2/rp_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/r_t - 1/a_f)) - sqrt(mu_s * (2/r_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (rp_b > r_t) then  ! If ra_b is larger, then it's better to do the inclination change post-apside change

        rstar = r_t / rp_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b))

        v2 = sqrt(dV1 ** 2 + dVi **2) + abs(dV2)

    else

        rstar = rp_b/r_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b))
        v2 = abs(dV1) + sqrt(dV2 ** 2 + dVi **2)

    end if

    cheap = min(v1, v2)

end subroutine cheapest


subroutine data_index(k, outval1, outval2)

    integer, intent(in) :: k
    integer, intent(out) :: outval1, outval2 

        if (k <= 2791440) then

        outval1 = floor(REAL(k)/360)*360
        outval2 = outval1 + 360

        else if (k > 2791440 .AND. k <= 2800674) then

        outval1 = floor((REAL(k)-2791440)/32)*32+2791440
        outval2 = outval1+32

        else if (k > 2800674) then

        outval1 = floor((REAL(k)-2800674)/360)*360+2800674
        outval2 = outval1+360

    end if

end subroutine data_index

!
! Newton-Raphson method to solve Kepler's equation
!

subroutine kepler_eq(x0, e, M, output, success)

    !
    ! Newton-Raphson method for the solution of Kepler's equation
    !
    ! Subroutine requires obj_fun_kep() function to proceed.
    !

    real, intent(in)     :: M           ! Value of M to solve for
    real, intent(in)     :: x0          ! Initial guess
    real, intent(in)     :: e           ! Eccentricity

    real, intent(out)    :: output      ! Result (eccentric anomaly)
    
    logical, intent(out) :: success     ! Success flag

    real, parameter      :: tol = 1e-09 ! Solution tolerance (may need tweaking)
    real, parameter      :: eps = 1e-06 ! Newton-Raphson tolerance

    real                 :: fx1, fx2    ! Function evaluations
    real                 :: fprime      ! Derivative estimation
    real                 :: x           ! Previous guess of E
    real                 :: xnew        ! New guess for E

    integer              :: i           ! Loop variable
    integer              :: maxiter     ! Maximum number of iterations: default = 50

    result = 0.0
    success = .false.
    x = x0

    maxiter = 50

    do i = 1, maxiter

        fx1 = obj_fun_kep(x, e, M)      ! Compute two guesses for the solution
        fx2 = obj_fun_kep(x+eps, e, M)

        fprime = ( fx2 - fx1 ) / eps    ! Estimate derivative

        xnew = x - fx1/fprime

        if ( abs(x-xnew) <= tol ) then

            success = .true.            ! Process was successful
            output = xnew               ! Leave the loop early
            exit

        end if


        x = xnew

    end do

end subroutine kepler_eq

!
! Newton-Raphson objective function
!

real function obj_fun_kep(x, e, M)

    real, intent(in)  :: x 
    real, intent(in)  :: M 
    real, intent(in)  :: e

    obj_fun_kep = x - e * sin(x) - M

end function

!
! Wrapper routine for converting mean anomaly to true anomaly
!

subroutine mean_to_true(M, e, nu_t)

    real, intent(in)  :: M              ! Mean anomaly
    real, intent(in)  :: e              ! Eccentricity
    
    real, intent(out) :: nu_t           ! True anomaly 

    real              :: E_t            ! Target eccentric anomaly

    logical           :: success        ! N-R iteration success flag

    ! Use kepler_eq() to get the equivalent eccentric anomaly 
    ! for this problem

    call kepler_eq(M, e, M, E_t, success)

    if (success .eqv. .false.) then

        error stop "Kepler's equation could not be solved."
    
    end if

    !
    ! Use intrinsic inverse tan to compute nu_t
    ! Uses atan2 to attempt to preserve quadrant
    !

    nu_t = 2. * atan2(sqrt((1+e)/(1-e)) * tan(E_t/2), 1.)

end subroutine mean_to_true
