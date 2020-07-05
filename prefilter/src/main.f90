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
!   METADATA
!   --------
!   Date: 31/03/18
!   Author: Jack Tyler
!
!   CHANGELOG
!   ---------
!   1/1/18: Bug fixes
!   2/2/2019: Minor refactoring
!   4/7/2020: Doxygen
!


program main

    use constants
    use problem_parameters
    use variable_init

    implicit none

    call mpi_variable_initialisation()

    call pre_filter()

    call variable_destruct()

end program main

subroutine pre_filter()

    use precision_kinds
    use constants
    use problem_parameters
    use variable_init 

    integer                                     :: i, k                                                  ! Loop variables 
    integer                                     :: bodyID                                                ! SPK ID of the body being investigated
    integer                                     :: week_id                                               ! Week ID to compute the solution for
                                                                                                         ! For use on Lyceum, the input was parameterised by date. At higher core counts,
                                                                                                         ! It becomes feasible to do all dates at once.
    integer                                     :: file_unit
    integer                                     :: week_index

    double precision                            :: rp_b                                                  ! SMA of body
    double precision                            :: a_b
    double precision                            :: e_b                                                   ! Eccentricity of body 
    double precision                            :: in_b                                                  ! Inclination of body (rad!!)
    double precision                            :: rp_t                                                  ! Periapse radius of target
    double precision                            :: e_t
    double precision                            :: in_t                                                  ! Inclination of target (rad)
    
    double precision                            :: mindV                                                 ! minimum delta V
    double precision                            :: epoch

    character(len=260)                          :: bodystring                                            ! String of the body SPK ID (file IO), OUTPUT_FILE for directory prefix
    character(len=8)                            :: fmt = '(I7.7)'                                        ! Formatting for result I/O
    character(len=4)                            :: week_string

    !
    ! Loop through all dates to consider - each node takes one date
    !

    do week_index = (mpi_id_world+1), size(dataset), mpi_world_size

        !
        ! Get week id from database; read into string for opening
        !

        week_id = dataset(week_index)
        write(week_string, '(I4)') week_id                                                                       ! Ephemeris seconds

        ! String version of week ID for determining file

        call week_to_time(week_id, epoch)

        ! Open the results file for this week

        file_unit = 150 + mpi_id_world
        open(unit=file_unit, file='../results/week_'//trim(adjustl(week_string)))

        !
        ! Loop through all objects in the MBD
        !

        do i = 1, num_candidates ! All bodies

            !
            ! Get the body ID from the database and write to a string
            !

            bodyID = asteroid_database(i)
            write(bodystring, '(I7)') bodyID

            !
            ! Load the candidate ephemeris; get the state at the desired time; convert to OEs; unload
            !

            call FURNSH(candidate_ephem_prefix//trim(adjustl(bodystring))//'.bsp')
            call SPKEZR(bodystring, epoch, 'ECLIPJ2000', 'NONE', 'Sun', asteroid_state, dum)
            call OSCELT(asteroid_state, epoch, mu, asteroid_elements)
            call UNLOAD(candidate_ephem_prefix//trim(adjustl(bodystring))//'.bsp')

            !
            ! Extract individual OEs from the array
            !

            rp_b = asteroid_elements(1)
            a_b  = rp_b / (1.d0 - e_b)                                                                            ! a is non-dim: convert to km
            e_b  = asteroid_elements(2)
            in_b = asteroid_elements(3)

            !
            ! Loop through all targets in the database
            !

            exitloop: do k = 1, num_orbits

                !
                ! Extract target orbital elements; dimensionalize
                !
                
                rp_t = manifold_data_elements(k,1)*au                                                              ! rp_t is non-dim: convert to km
                e_t  = manifold_data_elements(k,2)
                in_t = manifold_data_elements(k,3)

                !
                ! Compute the cheapest transfer to this condition
                !

                call cheapest(a_b, e_b, in_b, rp_t, e_t, in_t, mindV)                                               ! Compute if we have any possible transfers

                !
                ! If transfer less than the threshold, add to the output file
                !

                if (mindV < .7d0) then

                    call file_add(bodyID, file_unit, mindV)                                                         ! Open results file, add the SPK ID in, and close again
                    candidate_count = candidate_count + 1
                    exit exitloop ! Exit the target loop; we don't need to look any further

                end if

            end do exitloop

            if (mod(i, 1000) .eq. 0) then

                write(*,*) "Completed processing object number ", i
                write(*,*) "Candidates found: ", candidate_count

            end if

        end do

        close(unit=file_unit)

    end do

    close(file_unit)

end subroutine pre_filter

!
! Routine to add the SPK ID of a body to a file if it is deemed to be a candidate
!

subroutine file_add(bodyID, unitnum, mindv)

    use precision_kinds

    ! Subroutine intents

    integer,          intent(in)    :: bodyID
    integer,          intent(in)    :: unitnum
    integer,          intent(in)    :: k

    double precision, intent(in)    :: mindv

    write(unitnum, *) bodyID, mindv
    close(unitnum)

end subroutine file_add

!
! Gets the Hohmann transfer velocity for a given set of OEs
!

subroutine get_transfer_velocity(a_0, a_f, r_t, r_b, in_t, in_b, transfer_vel)

    use precision_kinds
    use constants
    
    implicit none
    
    double precision, intent(in)    :: a_0                              ! Initial SMA  (km)
    double precision, intent(in)    :: a_f                              ! Final SMA (km)
    double precision, intent(in)    :: r_t                              ! Target initial radius (km)
    double precision, intent(in)    :: r_b                              ! Candidate initial radius(km)
    double precision, intent(in)    :: in_t                             ! Inclination of the target (rad)
    double precision, intent(in)    :: in_b                             ! Inclination of the candidate (rad)
    
    double precision, intent(out)   :: transfer_vel                     ! Transfer velocity (km/s)
    
    double precision                :: rstar                            ! Quotient for inclination change calculation
    double precision                :: dVi                              ! Velocity required for inclination change (km/s)
    double precision                :: dV1                              ! Velocity required for first Hohmann maneouvre (km/s)
    double precision                :: dV2                              ! Velocity required for second Hohmann maneouvre (km/s)
    double precision                :: a_int                            ! SMA of the intermediate transfer arc (km)

    !
    ! First, compute the intermediate and final semimajor axes
    !

    a_int = (r_t + r_b) / 2.

    !
    ! Compute the velocity required to change the size of the apses
    !

    dV1 = sqrt( mu * (2/r_b - 1/a_int) ) - sqrt( mu * (2/r_b - 1/a_0) )
    dV2 = sqrt( mu * (2/r_t - 1/a_f) ) -   sqrt( mu * (2/r_t - 1/a_int) )

    !
    ! Now the inclination change: only do it at the point where it is cheapest.
    ! Below is taken from Sanchez et. al., 2016b.
    !

    if (r_t > r_b) then  ! If r_t is larger, then it's better to do the inclination change post-apside change

        rstar = r_b / r_t
        dVi = 2 * sqrt( (mu/a_int) * rstar ) * sin( abs( in_t - in_b ) / 2.0 )  ! Inclination change
        transfer_vel = abs(dV1) + sqrt(dVi ** 2 + dV2 ** 2)

    else

        rstar = r_t / r_b
        dVi = 2 * sqrt( (mu/a_int) * rstar) * sin( abs( in_t - in_b ) / 2.0 )   ! Inclination change
        transfer_vel = sqrt(dV1 ** 2 + dVi ** 2) + abs(dV2)

    end if

end subroutine get_transfer_velocity

!
! Obtain the lowest Hohmann transfer velocity for the given target and candidate positions
!

subroutine cheapest(a_b, e_b, in_b, rp_t, e_t, in_t, cheap)

    use precision_kinds
    use constants
    
    implicit none

    double precision, intent(in)    :: a_b                  ! SMA of the candidate (km)
    double precision, intent(in)    :: e_b                  ! Eccentricity of the candidate
    double precision, intent(in)    :: in_b                 ! Inclination of the candiate (rad)
    double precision, intent(in)    :: rp_t                 ! Radius of periapsis of the target (km)
    double precision, intent(in)    :: e_t                  ! Eccentricity of the target
    double precision, intent(in)    :: in_t                 ! Inclination of the target (rad)

    double precision, intent(out)   :: cheap                ! Cheapest Hohmann transfer (km/s)

    double precision                :: ra_b                 ! Radius of apoapsis of candidate (km)
    double precision                :: rp_b                 ! Radius of periapsis of candidate (km)
    double precision                :: ra_t                 ! Radius of apoaosis of target (km)
    double precision                :: a_0                  ! SMA of candidate (km)
    double precision                :: a_t                  ! SMA of target (km)
    double precision                :: a_f                  ! Equivalent SMA of target (doesn't really exist)
    double precision                :: v1                   ! Total velocity required for first maneouvre (km/s)
    double precision                :: v2                   ! Total velocity required for second maneouvre (km/s)
    double precision                :: v3                   ! Total velocity required for first maneouvre (km/s)
    double precision                :: v4                   ! Total velocity required for second maneouvre (km/s)

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

!
! Convert an integer ID week to ephemeris seconds (assumes baseline is Jan 1, 2025)
!

subroutine week_to_time(week, time)

    integer, intent(in)             :: week                 ! Integer ID of week being studied

    double precision, intent(out)   :: time                 ! Ephemeris time of week being studied

    double precision                :: baseline_time

    call FURNSH('../data/naif0008.tls')

    call STR2ET('Jan 1, 2025 00:00', baseline_time)         ! Baseline time

    time = baseline_time + week * 86400.d0 * 7.d0

end subroutine week_to_time
