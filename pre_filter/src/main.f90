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
    use constants                                                                                        ! Standardises constants for use in all Fortran codes
    use problem_parameters
    use variable_init

    implicit none

    call variable_initialisation()

    call pre_filter()

    call variable_destruct()

end program main

subroutine pre_filter()

    use precision
    use constants
    use problem_parameters
    use variable_init 

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
    real(kind=dp)                               :: in_t                                                  ! Inclination of target (rad)
    
    real(kind=dp)                               :: mindV                                                 ! minimum delta V
    real(kind=dp)                               :: t1, t2, t3, t4, t5, date, o, Om, M, p                 ! Temporary variables
    real(kind=dp)                               :: ephemeris_time

    character(len=260)                          :: bodystring                                            ! String of the body SPK ID (file IO), OUTPUT_FILE for directory prefix
    character(len=8)                            :: fmt = '(I7.7)'                                        ! Formatting for result I/O

    ! Variable initialisations

    ! Generate the orbital elements of the candidate at the desired time

    call getarg(1, desireddate)                                                                          ! Ephemeris seconds

    do i = 1, num_targets ! All bodies

        ! Get the body ID from the database and write to a string

        bodyID = database(i)
        write(bodystring, *) bodyID

        call FURNSH(candidate_ephem_prefix//bodyID//'.bsp')
        call SPKEZR(bodyID, desireddate, 'ECLIPJ2000', 'NONE', 'Sun', asteroid_state, dum)
        call OSCELT(asteroid_state, desireddate, mu, asteroid_elements)

        a_b  = asteroid_elements(1)                                                                            ! a is non-dim: convert to km
        e_b  = asteroid_elements(2)
        in_b = asteroid_elements(3)

        exitloop: do k = 1, num_targets
            
            rp_t = manifold_data_elements(k,1)*au                                                               ! rp_t is non-dim: convert to km
            e_t  = manifold_data_elements(k,2)
            in_t = manifold_data_elements(k,3)

            call cheapest(a_b, e_b, in_b, rp_t, e_t, in_t, mindV)                                               ! Compute if we have any possible transfers

            if (mindV < .7) then 

                call file_add(bodyID, unitno, mindV, k)                                                         ! Open results file, add the SPK ID in, and close again
                candidate_count = candidate_count + 1
                exit exitloop                                                                                   ! If we have a candidate, then we don't care about any other transfer opportunities; stop.

            end if

        end do

        if (mod(i, 10000) .eq. 0) then

            write(*,*) "Completed processing object number ", i
            write(*,*) "Candidates found: ", candidate_count

        end if

        call UNLOAD(candidate_ephem_prefix//bodyID//'.bsp')

    end do

    write(*,*) "==================== Scan completed ========================"
    write(*,*) "Objects scanned: ", itercount
    write(*,*) "Candidates found: ", candidate_count
    write(*,*) "Compute time: ", finish-start, "seconds."

end subroutine pre_filter

!
! Routine to add the SPK ID of a body to a file if it is deemed to be a candidate
! WARNING: May have issues if multiple threads attempt to I/O at once.
!

subroutine file_add(bodyID, unitnum, mindv, k)

    ! Subroutine intents

    integer,        intent(in) :: bodyID
    integer,        intent(in) :: unitnum
    integer,        intent(in) :: k

    real(kind=dp),  intent(in) :: mindv

    logical                     :: exist

    inquire(file=OUTPUT_FILE, exist=exist)

    if (exist) then

      open(unitnum, file=OUTPUT_FILE, status="old", position="append", action="write")

    else

      open(unitnum, file=OUTPUT_FILE, status="new", action="write")

    end if

    ! Write the body, the mindv, and the family (bounds)

    write(unitnum, *) bodyID, mindv
    close(unitnum)

end subroutine file_add

!
! Subroutine to check if the results file already exists in the working folder
!

subroutine file_check()

    logical           :: exist

    OUTPUT_FILE = 'results.txt'

    inquire(file=OUTPUT_FILE, exist=exist)

    if (exist) then

        open(20, file=OUTPUT_FILE)
        close(20, status='delete')

    end if

end subroutine file_check


subroutine get_transfer_velocity(a_0, a_f, r_t, r_b, in_t, in_b, transfer_vel)

    use precision
    use constants
    
    implicit none
    
    real(kind=dp), intent(in)       :: a_0                              ! Initial SMA  (km)
    real(kind=dp), intent(in)       :: a_f                              ! Final SMA (km)
    real(kind=dp), intent(in)       :: r_t                              ! Target initial radius (km)
    real(kind=dp), intent(in)       :: r_b                              ! Candidate initial radius(km)
    real(kind=dp), intent(in)       :: in_t                             ! Inclination of the target (rad)
    real(kind=dp), intent(in)       :: in_b                             ! Inclination of the candidate (rad)
    
    real(kind=dp), intent(out)      :: transfer_vel                     ! Transfer velocity (km/s)
    
    real(kind=dp)                   :: rstar                            ! Quotient for inclination change calculation
    real(kind=dp)                   :: dVi                              ! Velocity required for inclination change (km/s)
    real(kind=dp)                   :: dV1                              ! Velocity required for first Hohmann maneouvre (km/s)
    real(kind=dp)                   :: dV2                              ! Velocity required for second Hohmann maneouvre (km/s)

    
    ! First, compute the intermediate and final semimajor axes

    a_int = (r_t + r_b) / 2.

    ! Compute the velocity required to change the size of the apses

    dV1 = sqrt( mu * (2/r_b - 1/a_int) ) - sqrt( mu * (2/r_b - 1/a_b) )
    dV2 = sqrt( mu * (2/r_t - 1/a_f) ) -   sqrt( mu * (2/r_t - 1/a_int) )

    !
    ! Now the inclination change: only do it at the point where it is cheapest
    !

    if (r_t > r_b) then  ! If r_t is larger, then it's better to do the inclination change post-apside change

        rstar = r_b / r_t
        dVi = 2 * sqrt( (mu/a_int) * rstar ) * sin( abs( in_t - in_b ) / 2.0 )
        transfer_vel = abs(dV1) + sqrt(dVi ** 2 + dV2 ** 2)

    else

        rstar = r_t / r_b
        dVi = 2 * sqrt( (mu/a_int) * rstar) * sin( abs(in_t - in_b) / 2.0 )
        transfer_vel = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

end subroutine get_transfer_velocity


subroutine cheapest(a_b, e_b, in_b, rp_t, e_t, in_t, M_t, cheap)

    implicit none

    use precision
    use constants

    real(kind=dp), intent(in)       :: a_b                  ! SMA of the candidate (km)
    real(kind=dp), intent(in)       :: e_b                  ! Eccentricity of the candidate
    real(kind=dp), intent(in)       :: in_b                 ! Inclination of the candiate (rad)
    real(kind=dp), intent(in)       :: rp_t                 ! Radius of periapsis of the target (km)
    real(kind=dp), intent(in)       :: e_t                  ! Eccentricity of the target
    real(kind=dp), intent(in)       :: in_t                 ! Inclination of the target (rad)
    real(kind=dp), intent(in)       :: M_t                  ! Mean anomaly of the target (rad)

    real(kind=dp), intent(out)      :: cheap                ! Cheapest Hohmann transfer (km/s)

    real(kind=dp)                   :: ra_b                 ! Radius of apoapsis of candidate (km)
    real(kind=dp)                   :: rp_b                 ! Radius of periapsis of candidate (km)
    real(kind=dp)                   :: ra_t                 ! Radius of apoaosis of target (km)
    real(kind=dp)                   :: a_t                  ! SMA of target (km)
    real(kind=dp)                   :: a_int                ! SMA of intermediate transfer arc (km)
    real(kind=dp)                   :: a_f                  ! Equivalent SMA of target (doesn't really exist)
    real(kind=dp)                   :: v1                   ! Total velocity required for first maneouvre (km/s)
    real(kind=dp)                   :: v2                   ! Total velocity required for second maneouvre (km/s)
    real(kind=dp)                   :: v3                   ! Total velocity required for first maneouvre (km/s)
    real(kind=dp)                   :: v4                   ! Total velocity required for second maneouvre (km/s)

    ! First, Compute apoapsis and periapsis distance for the body
    
    ra_b = a_b * (1 + e_b)
    rp_b = a_b * (1 - e_b)

    a_t = rp_t / (1 - e_t)
    ra_t = a_t * (1 + e_t)

    a_0 = a_b
    a_f = a_t

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
    ! Case 1, ap -> ap (ra_b + ra_t)
    !

    call get_transfer_velocity(a_0, a_f, ra_t, ra_b, in_t, in_b, v1)
    
    !
    ! Case 2, ap -> peri (ra_b + rp_t)
    !

    call get_transfer_velocity(a_0, a_f, ra_b, rp_t, in_t, in_b, v2)

    !
    ! Case 3, peri -> ap (rp_b + ra_t)
    !

    call get_transfer_velocity(a_0, a_f, rp_b, ra_t, in_t, in_b, v3)

    !
    ! Case 4, peri -> peri (rp_b + rp_t)
    !

    call get_transfer_velocity(a_0, a_f, rp_b, rp_t, in_t, in_b, v4)
    
    cheap = min(v1, v2, v3, v4)

end subroutine cheapest
