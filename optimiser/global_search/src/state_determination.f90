! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     CANDIDATE STATE DETERMINATION
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The following code implements the technique
! outlined by Sanchez. et. al. to determine the
! state of a candidate based on its
! distance from the Earth about a full synodic
! period
!
! de414.bsp, naif0008.tls, and the candidate ephemeris
! must be present in ../data.
!

! @brief Contains routines for obtaining and manipulating the states of NEOs
module state_determination

    implicit none

        contains

            ! @brief DEPRECATED: Implements method following Sanchez et. al, 2016 - gets the nominal state of the NEO in first synodic period NEO outside of pi/8 cone
            ! @param[in] targ_can SPK ID of asteroid to find nominal state for
            ! @param[out] time_lower Beginning of synodic period
            ! @param[out] time_upper End of synodic period
            subroutine get_state(targ_can, time_lower, time_upper)
                use constants

                implicit none

                CHARACTER(*),     intent(in)    :: targ_can                                 ! String of the candidate SPK ID

                double precision, intent(out)   :: time_lower                                ! Lower bound of the solution search time
                double precision, intent(out)   :: time_upper                                ! Upper bound of the solution search time

                double precision                :: epoch                                    ! Epoch: start time
                double precision                :: ang_can                                  ! Angle of the candidate with respect to the Sun
                double precision                :: tau_can                                  ! Period of the candidate
                double precision                :: epoch_upper                              ! Upper bound for search space
                double precision                :: last_checked                             ! Time candidate was last in pi/8
                double precision                :: dum                                      ! Dummy variable
                double precision                :: max_distance                             ! Maximum distance between the candidate and the Earth
                double precision, dimension(6)  :: state_can, state_syn                     ! States of the can.,Earth, and the synodic state of the candidate respectively

                CHARACTER(len=6)                :: abcorr                                   ! Abberation correcton string: initialised to NONE
                CHARACTER(len=5)                :: obs                                      ! Observing body string
                CHARACTER(len=12)               :: coord                                    ! Co-ordinate system of reference
                CHARACTER(len=19)               :: epoch_str                                ! String of the lower epoch bound
                CHARACTER(len=19)               :: epoch_upper_str                          ! String of the upper epoch bound
                CHARACTER(len=7)                :: targ_ear                                 ! String of the target body (Earth)

                integer                         :: istate = 0                               ! Success flag, iteration counter, max. array locator
                integer                         :: epoch_counter                            ! Epoch loop variable

                max_distance = 0                                                            ! Initialise to trivially low to force update
                
                ! Ephemeris options
                
                targ_ear        = 'Earth'
                abcorr          = 'NONE'
                obs             = 'Sun'
                coord           = 'ECLIPJ2000'
                epoch_str       = 'Jan 1, 2020 00:00'       
                epoch_upper_str = 'Jan 1, 2100 00:00'

                ! Get the lower bound of the epoch in terms of ephemeris seconds (et)

                call STR2ET(epoch_str, epoch)

                ! As above, but with the upper bound

                call STR2ET(epoch_upper_str, epoch_upper)

                ! Initialise 'trip' counter to be the start of the domain

                last_checked = epoch                                                    ! To prevent errors on first run

                ! Initialise the period to which we enforce outside of pi/8 (seconds)

                tau_can = 365.25d0 * 86400.d0                                           ! One synodic periodic, seconds

                do epoch_counter = int(epoch), int(epoch_upper), 43200                  ! Deprecated feature: non-integer loop sentinels

                    ! Get the initial state of the candidate

                    call SPKEZR(targ_can, epoch_counter, 'ECLIPJ2000', abcorr,&
                                obs, state_can, dum)

                    ! Rotate into synodic frame to check the angle w.r.t. Earth;
                    ! overwrite current state_can with new synodic version

                    call SYNODIC_ROTATE(state_can, epoch_counter*1.d0, state_syn)

                    ! Compute angle of candidate w.r.t. Earth

                    ang_can = datan2(state_syn(2), state_syn(1))

                    ! If we are within the pi/8 plane, reset the trip counter

                    if (abs(ang_can) < pi/8) then

                        last_checked = epoch_counter

                    end if

                    ! If we have been outside of the pi/8 plane for more than tau_can, exit
                    ! epoch_counter is always >= last_checked

                    if ((epoch_counter-last_checked) > tau_can) then

                        istate=1                                                    ! Switch flag
                        exit                                                        ! Exit the loop

                    end if

                end do

                ! Check if we found a solution

                if (istate /= 1) then

                    error stop "Solution not found. Stopping..."

                end if

                time_lower = last_checked
                time_upper = last_checked+tau_can

            end subroutine get_state

            ! @brief Wrapper function to convert a calendar date string to ephemeris seconds
            ! @param[in] str_in Calendar date string to adjust
            ! @param[out] time_out Equivalent time in ephemeris seconds
            subroutine time_str(str_in, time_out)

                CHARACTER(*), intent(in)        :: str_in

                double precision, intent(out)   :: time_out
                
                call STR2ET(str_in, time_out)

            end subroutine time_str

            ! @brief Get the position of a candidate NEO from the ephemeris
            ! @param[in] t_epoch Time to obtain state at in ephemeris seconds
            ! @param[out] state_can State of the candidate at t_epoch
            ! Candidate SPK ID is taken from the problem parameters module.
            subroutine CANDIDATE_POSITION(t_epoch, state_can)
                 
                use problem_parameters                                                  ! Gives access to targ_can

                double precision, intent(in)   :: t_epoch                               ! Initial time
                double precision, intent(out)  :: state_can(6)                          ! Final state
                
                CHARACTER(len=6)               :: abcorr                                ! Abberation correcton string: initialised to NONE
                CHARACTER(len=5)               :: obs                                   ! Observing body string
                CHARACTER(len=12)              :: coord    
                
                double precision               :: dum

                abcorr = 'NONE'
                obs    = 'Sun'
                coord  = 'ECLIPJ2000'

                ! Get position from ephemeris
                call SPKEZR(targ_can, t_epoch, 'ECLIPJ2000', 'NONE', &
                            'Sun', state_can, dum)

            end subroutine CANDIDATE_POSITION

            ! @brief Rotate a state from the inertial frame into the synodic frame and non-dimensionalise the resulting state
            ! @param[in] state_in State to rotate; epoch of J2000 in ephemeris seconds, km/s
            ! @param[in] t Time at which to rotate; ephemeris seconds
            ! @param[out] state_out Rotated state; dimensionless
            subroutine SYNODIC_ROTATE(state_in, t, state_out)
                 
                use constants

                double precision, intent(in)   :: state_in(6)
                double precision, intent(in)   :: t

                double precision, intent(out)  :: state_out(6)

                double precision, parameter    :: theta_0 = 100.3762 * pi / 180.d0          ! Angle of the Earth at J2000

                double precision               :: r(3), rdot(3)
                double precision               :: t_ir(3, 3)                                ! Rotation matrix
                double precision               :: t_irdot(3, 3)                             ! Derivative of rotation matrix
                double precision               :: total_angle                               ! Total angle to rotate
                double precision               :: cang                                      ! Cosine of the angle
                double precision               :: sang                                      ! Sine of the angle
                
                ! Get total angle - OK to make large as cos/sin functions will modulo 2pi the answer
                total_angle = theta_0 + (t * time_to_angle_quotient)

                ! Pre-compute cos/sin of angle
                cang = cos(total_angle)
                sang = sin(total_angle)

                ! Non-dimensionalise first
                r = state_in(1:3) * position_non_dimensionalise_quotient
                rdot = state_in(4:6) * velocity_non_dimensionalise_quotient

                ! Construct T
                t_ir(1,1) = cang; t_ir(1, 2) = sang; t_ir(1,3) = 0.0
                t_ir(2,1) = -sang; t_ir(2,2) = cang; t_ir(2,3) = 0.0
                t_ir(3,1) = 0.0; t_ir(3,2) = 0.0; t_ir(3,3) = 1.0

                ! Construct T_dot
                t_irdot(1,1) = -sang; t_irdot(1,2) = cang; t_irdot(1,3) = 0.0
                t_irdot(2,1) = -cang; t_irdot(2,2) = -sang; t_irdot(2,3) = 0.0
                t_irdot(3,1) = 0.0; t_irdot(3,2) = 0.0; t_irdot(3,3) = 0.0

                ! Multiply
                state_out(1:3) = matmul(t_ir, r)
                state_out(4:6) = matmul(t_ir, rdot) + matmul(t_irdot, r)

                ! Shift barycenter
                state_out(1) = state_out(1) - mu3bp

            end subroutine SYNODIC_ROTATE

            ! @brief Rotate a state from the synodic frame (dimensionless) into the inertial frame
            ! @param[in] state_in State to rotate; epoch of J2000 in ephemeris seconds, dimensionless
            ! @param[in] t Time at which to rotate; ephemeris seconds
            ! @param[out] state_out Rotated state; dimensionless
            subroutine GLOBAL_ROTATE(state_in, t, state_out)

                use constants

                double precision, intent(in)    :: state_in(6)
                double precision, intent(in)    :: t

                double precision, intent(out)   :: state_out(6)

                double precision, parameter     :: theta_0 = 100.3762 * 3.141592 / 180.d0           ! Angle of the Earth at J2000

                double precision                :: r(3), rdot(3)
                double precision                :: t_ir(3, 3)
                double precision                :: t_irdot(3, 3)
                double precision                :: state_translate(6)
                double precision                :: total_angle
                double precision                :: cang
                double precision                :: sang 
                
                ! Get total angle - OK to make large as cos/sin functions will modulo 2pi the answer
                total_angle = theta_0 + t * time_to_angle_quotient

                ! Shift barycenter
                state_translate = state_in
                state_translate(1) = state_in(1) + mu3bp

                ! Pre-compute cos/sin of angle
                cang = cos(total_angle)
                sang = sin(total_angle)

                ! Non-dimensionalise first
                r = state_translate(1:3)
                rdot = state_translate(4:6)

                ! Construct T
                t_ir(1,1) = cang; t_ir(1, 2) = sang; t_ir(1,3) = 0.0
                t_ir(2,1) = -sang; t_ir(2,2) = cang; t_ir(2,3) = 0.0
                t_ir(3,1) = 0.0; t_ir(3,2) = 0.0; t_ir(3,3) = 1.0

                ! Construct T_dot
                t_irdot(1,1) = -sang; t_irdot(1,2) = cang; t_irdot(1,3) = 0.0
                t_irdot(2,1) = -cang; t_irdot(2,2) = -sang; t_irdot(2,3) = 0.0
                t_irdot(3,1) = 0.0; t_irdot(3,2) = 0.0; t_irdot(3,3) = 0.0

                ! Multiply
                t_ir = transpose(t_ir)
                t_irdot = transpose(t_irdot)

                state_out(1:3) = matmul(t_ir, r)
                state_out(4:6) = matmul(t_ir, rdot) + matmul(t_irdot, r)

            end subroutine GLOBAL_ROTATE

end module