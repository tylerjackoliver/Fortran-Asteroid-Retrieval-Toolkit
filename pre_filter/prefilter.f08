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
    use constants                                                                           ! Standardises constants for use in all Fortran codes
    use iso_fortran_env, only : error_unit 

    implicit none

    integer                                     :: nthreads = 2                                          ! Number of threads to use
    integer                                     :: candidatecount
    integer                                     :: itercount
    integer                                     :: unitno                                                ! Counting variables
    integer                                     :: i, j, k                                               ! Loop variables 
    integer                                     :: bodyID                                                ! SPK ID of the body being investigated
    integer                                     :: desireddate                                           ! Week ID to compute the solution for
                                                                                                         ! For use on Lyceum, the input was parameterised by date. At higher core counts,
                                                                                                         ! It becomes feasible to do all dates at once.

    real(kind=dp)                               :: a_b                                                   ! SMA of body
    real(kind=dp)                               :: e_b                                                   ! Eccentricity of body 
    real(kind=dp)                               :: in_b                                                  ! Inclination of body (rad!!)

    real(kind=dp)                               :: rp_t                                                  ! Periapse radius of target
    real(kind=dp)                               :: e_t
    real(kind=dp)                               :: in_t                                                  ! Inclination of target (deg!!)

    real(kind=dp)                               :: start                                                 ! OMP WTIME timing
    real(kind=dp)                               :: finish                                                ! ""
    
    real(kind=dp)                               :: M_t                                                   ! Mean anomaly of the target
    real(kind=dp)                               :: mindV                                                 ! minimum delta V
    real(kind=dp)                               :: t1, t2, t3, t4, t5, date, o, Om, M, p                 ! Temporary variables

    integer, dimension(17719)                   :: database                                              ! SPK ID Database

    real(kind=dp),    dimension(1353, 3)        :: bodyData                                              ! rp_b, e_b, inb,  of the target body (needed only for Hohmann)
    real(kind=dp),    dimension(360000,4)       :: targetData                                            ! rp_b, e_t, in_t, M_t of the target manifold section

    character(len=260)                          :: bodystring, filepath                                  ! String of the body SPK ID (file IO), filepath for directory prefix
    character(len=8)                            :: fmt = '(I7.7)'                                        ! Formatting for result I/O

    ! Variable initialisations

    candidatecount = 0                                                                        ! Number of candidates found                                                          
    itercount = 0                                                                             ! Number of weeks processed
    filepath = 'corrected_oe/'                                                                ! Directory prefix: alter as req'd

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

    call OMP_SET_NUM_THREADS(nthreads)    

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

    open(unit=99,file="successful_OEs.txt")

    ! Write progress message and begin timing

    write(*,*) "Loading complete. Beginning compute sequence for body #", desireddate
    
    start = omp_get_wtime()

    !$OMP PARALLEL DO PRIVATE(I,J,K, bodyid, bodystring, unitno, mindV, bodyData, a_b, e_b, inb, rp_t, e_t, in_t)
    do i = 1, size(database) ! All bodies

        ! Get the body ID from the database and write to a string

        bodyID = database(i)

        ! Unique unit number for each thread

        unitno = OMP_GET_THREAD_NUM() + 30

        ! Open file, read the input data and close again

        open(unit=unitno, file=trim(filepath)//trim(bodystring)//'_OE_corr.csv')
    
        do j = 1,1353

            read(unitno,*) date, bodyData(j,1), bodyData(j,2), bodyData(j,3), o, O, M, p

        end do

        close(unitno)

        ! Loop at the desired data through every body, and every target, to see if we have any possible transfers

        exitloop:  do j = 221,1353

            a_b = bodyData(j,1)*au                                                                                  ! a is non-dim: convert to km
            e_b = bodyData(j,2)
            inb = bodyData(j,3)                                                                                     ! Note: degrees - is converted in cheapest subroutine

            do k = 1, 360000
                
                rp_t = targetData(k,1)*au                                                                           ! rp_t is non-dim: convert to km
                e_t  = targetData(k,2)
                in_t = targetData(k,3)
                M_t  = targetData(k,4)

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

    finish = omp_get_wtime()

    write(*,*) "==================== Scan completed ========================"
    write(*,*) "Objects scanned: ", itercount
    write(*,*) "Candidates found: ", candidatecount
    write(*,*) "Compute time: ", finish-start, "seconds."
    close(99)

end program main

!
! Routine to add the SPK ID of a body to a file if it is deemed to be a candidate
! WARNING: May have issues if multiple threads attempt to I/O at once.
!

subroutine file_add(bodyID, unitnum, mindv, k)

    ! Subroutine intents

    character(len=13)   :: filepath ! Define address of file to be read

    integer, intent(in) :: bodyID
    integer, intent(in) :: unitnum
    integer, intent(in) :: k

    real(kind=dp), intent(in)    :: mindv
    
    integer             :: outval1, outval2

    logical             :: exist

    ! call where_from(k, outstring)

    filepath = "results.txt"

    inquire(file=filepath, exist=exist)

    if (exist) then

      open(unitnum, file=filepath, status="old", position="append", action="write")

    else

      open(unitnum, file=filepath, status="new", action="write")

    end if

    ! Write the body, the mindv, and the family (bounds)

    write(unitnum, *) bodyID, mindv, outval1, outval2
    close(unitnum)

end subroutine file_add

!
! Subroutine to check if the results file already exists in the working folder
!

subroutine file_check()

    character(len=13) :: filepath
    logical           :: exist

    filepath = 'results.txt'

    inquire(file=filepath, exist=exist)

    if (exist) then

        open(20, file=filepath)
        close(20, status='delete')

    end if

end subroutine file_check


subroutine cheapest(a_b, e_b, inb, rp_t, e_t, in_t, M_t, cheap)

    real(kind=dp), intent(in)    :: a_b                  ! SMA of the candidate (km)
    real(kind=dp), intent(in)    :: e_b                  ! Eccentricity of the candidate
    real(kind=dp), intent(in)    :: inb                  ! Inclination of the candiate (deg)
    real(kind=dp), intent(in)    :: rp_t                 ! Radius of periapsis of the target (km)
    real(kind=dp), intent(in)    :: e_t                  ! Eccentricity of the target
    real(kind=dp), intent(in)    :: in_t                 ! Inclination of the target (rad)
    real(kind=dp), intent(in)    :: M_t                  ! Mean anomaly of the target (rad)

    real(kind=dp), intent(out)   :: cheap                ! Cheapest Hohmann transfer (km/s)

    real(kind=dp)                :: mu_s                 ! Solar graviational parameter (km ^ 2/s ^ 3)
    real(kind=dp)                :: ra_b                 ! Radius of apoapsis of candidate (km)
    real(kind=dp)                :: rp_b                 ! Radius of periapsis of candidate (km)
    real(kind=dp)                :: ra_t                 ! Radius of apoaosis of target (km)
    real(kind=dp)                :: a_t                  ! SMA of target (km)
    real(kind=dp)                :: a_int                ! SMA of intermediate transfer arc (km)
    real(kind=dp)                :: a_f                  ! Equivalent SMA of target (doesn't real(kind=dp)ly exist)
    real(kind=dp)                :: rstar                ! Quotient for inclination change calculation
    real(kind=dp)                :: dVi                  ! Velocity required for inclination change (km/s)
    real(kind=dp)                :: dV1                  ! Velocity required for first Hohmann maneouvre (km/s)
    real(kind=dp)                :: dV2                  ! Velocity required for second Hohmann maneouvre (km/s)
    real(kind=dp)                :: v1                   ! Total velocity required for first maneouvre (km/s)
    real(kind=dp)                :: v2                   ! Total velocity required for second maneouvre (km/s)
    real(kind=dp)                :: v3                   ! Total velocity required for first maneouvre (km/s)
    real(kind=dp)                :: v4                   ! Total velocity required for second maneouvre (km/s)

    double precision    :: pi = 4.*datan(1.d0)

    mu_s = 1.327124400189e11

    ! Convert inb to radian

    in_b = inb*pi/180.0

    ! First, Compute apoapsis and periapsis distance for the body
    
    ra_b = a_b * (1 + e_b)
    rp_b = a_b * (1 - e_b)

    a_t = rp_t / (1 - e_t)
    ra_t = a_t * (1 + e_t)

    !
    ! There are four possible cases of transfers:
    !       1. Ap of body -> target apoapsis
    !       2. Ap of body -> target periapsis
    !       3. Peri of body -> target apoapsis
    !       4. Peri of body -> target periapsis
    !
    ! We can, however, automatically eliminate two of these cases by simply choosing the most efficient point
    ! to do the inclination change. This is at the point at which the velocity of the object is *lowest*
    !

    !
    ! Case 1, ap -> ap
    !

    ! First, compute the intermediate and final semimajor axes

    a_0 = a_b
    a_int = (ra_b + ra_t) / 2.
    a_f = a_t

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/ra_b - 1/a_int)) - sqrt(mu_s * (2/ra_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/ra_t - 1/a_f)) - sqrt(mu_s * (2/ra_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (ra_t > ra_b) then  ! If ra_t is larger, then it's better to do the inclination change post-apside change

        rstar = ra_b / ra_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.0)
        v1 = abs(dV1) + sqrt(dVi ** 2 + dV2 ** 2)

    else

        rstar = ra_t/ra_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.0)
        v1 = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

    !
    ! Case 2, ap -> peri
    !

    ! Update the intermediate SMA

    a_int = (ra_b + rp_t) / 2.

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/ra_b - 1/a_int)) - sqrt(mu_s * (2/ra_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/rp_t - 1/a_f)) - sqrt(mu_s * (2/rp_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (rp_t > ra_b) then  ! If ra_b is larger, then it's better to do the inclination change post-apside change

        rstar = ra_b / rp_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.)
        v2 = abs(dV1) + sqrt(dV2 ** 2 + dVi **2)

    else

        rstar = rp_t/ra_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.)
        v2 = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

    !
    ! Case 3, peri -> ap
    !

    ! First, compute the intermediate and final semimajor axes

    a_0 = a_b
    a_int = (rp_b + ra_t) / 2.
    a_f = a_t

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/rp_b - 1/a_int)) - sqrt(mu_s * (2/rp_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/ra_t - 1/a_f)) - sqrt(mu_s * (2/ra_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (ra_t > rp_b) then  ! If ra_t is larger, then it's better to do the inclination change post-apside change

        rstar = rp_b / ra_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.0)
        v3 = abs(dV1) + sqrt(dV2 ** 2 + dVi ** 2)

    else

        rstar = ra_t/rp_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.0)
        v3 = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

    !
    ! Case 4, peri -> peri
    !

    ! Update the intermediate SMA

    a_int = (rp_b + rp_t) / 2.

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt(mu_s * (2/rp_b - 1/a_int)) - sqrt(mu_s * (2/rp_b - 1/a_b))
    dV2 = sqrt(mu_s * (2/rp_t - 1/a_f)) - sqrt(mu_s * (2/rp_t - 1/a_int))

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (rp_t > rp_b) then  ! If ra_b is larger, then it's better to do the inclination change post-apside change

        rstar = rp_b / rp_t
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.)
        v4 = abs(dV1) + sqrt(dV2 ** 2 + dVi ** 2)

    else

        rstar = rp_t / rp_b
        dVi = 2 * sqrt((mu_s/a_int) * rstar) * sin(abs(in_t - in_b)/2.)
        v4 = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

    cheap = min(v1, v2, v3, v4)

end subroutine cheapest
