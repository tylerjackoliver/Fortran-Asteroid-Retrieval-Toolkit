!------------------------------------------------------------------------------
! Fortran Asteroid Retrieval Tool (FART) v1.0: state_determination module
!------------------------------------------------------------------------------
!
! MODULE: State determination
!
!> @author
!> Jack Tyler, University of Southampton
!
! DESCRIPTION: 
!> Contains routines to handle constraints and manpulate and transform state
!! vectors
!
! REVISION HISTORY:
! 01 Mar 2018 - Initial Version
! 01 Jul 2020 - Refactoring; add Doxygen support
!------------------------------------------------------------------------------


module state_determination

    implicit none

    public  :: get_state, synodic_rotate, global_rotate, candidate_position, rotator
    private :: time_str 

    contains

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Obtains the lower and upper values of time that constrain the transfer
        !! epoch as per Sanchez et. al., 2016. As of 2019, this routine is no longer used.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        ! 10 Oct 2019 - Deprecated
        !
        !> @param[in] targ_can
        !> @param[out] time_lower, time_upper
        !--------------------------------------------------------------------------- 

        subroutine get_state(targ_can, time_lower, time_upper)

            use constants

            implicit none

            double precision                :: epoch                                    ! Epoch: start time
            double precision                :: epoch_counter                            ! Epoch loop variable
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

            INTEGER                         :: istate = 0                               ! Success flag, iteration counter, max. array locator
            integer                         :: mpi_id_world
            integer                         :: mpi_err
            CHARACTER(*),  intent(in)       :: targ_can                                 ! String of the candidate SPK ID

            double precision, intent(out)  :: time_lower                                ! Lower bound of the solution search time
            double precision, intent(out)  :: time_upper                                ! Upper bound of the solution search time


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

            do epoch_counter = epoch, epoch_upper, 43200                            ! Deprecated feature: non-integer loop sentinels

                ! Get the initial state of the candidate

                call SPKEZR(targ_can, epoch_counter, 'ECLIPJ2000', abcorr,&
                            obs, state_can, dum)

                ! Rotate into synodic frame to check the angle w.r.t. Earth;
                ! overwrite current state_can with new synodic version

                call SYNODIC_ROTATE(state_can, epoch_counter, state_syn)

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

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Converts character array containing calendar date to ephemeris seconds.
        !! As of 2019, this routine is no longer used.
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        ! 10 Oct 2019 - Deprecated
        !
        !> @param[in] str_in
        !> @param[out] time_out
        !--------------------------------------------------------------------------- 

        subroutine time_str(str_in, time_out)

            CHARACTER(*), intent(in) :: str_in

            double precision, intent(out) :: time_out
            
            call STR2ET(str_in, time_out)

        end subroutine time_str

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Obtains the state of a candidate at a specified ephemeris time.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        !
        !> @param[in] t_epoch
        !> @param[out] state_can
        !--------------------------------------------------------------------------- 

        subroutine candidate_position(t_epoch, state_can)

            use problem_parameters

            double precision, intent(in)   :: t_epoch                               ! Initial time
            double precision, intent(out)  :: state_can(6)                          ! Final state

            double precision               :: dum

            call SPKEZR(targ_can, t_epoch, 'ECLIPJ2000', 'NONE', &
                        'Sun', state_can, dum)

        end subroutine candidate_position

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Rotates a given state in the J2000 frame about the +Z axis for an angle
        !! equivalent to the ephemeris time (e.g. an angle of 2*pi is 365.26*86400
        !! seconds.)
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        !
        !> @param[in] state_in, epoch, tt
        !> @param[out] state_out
        !--------------------------------------------------------------------------- 

        subroutine rotator(state_in, epoch, tt, state_out)
                
            use constants

            double precision, intent(in)   :: state_in(6)               ! Input state; to be rotated
            double precision, intent(in)   :: epoch                     ! Initial time; ephemeris seconds (seconds past J2000)
            double precision, intent(in)   :: tt                        ! Transfer time; set equal to zero to rotate only on epoch
            double precision, intent(out)  :: state_out(6)              ! Output state; rotated

            double precision               :: cang                      ! Pre-compute cosine of angle
            double precision               :: sang                      ! Pre-compute sine of angle

            double precision               :: theta                     ! Angle
            double precision               :: t_ir(3, 3)                ! Transformation matrix
            double precision               :: t_ir_dot(3, 3)            ! Derivative of the transformation matrix

            theta = (tt + epoch) * time_to_angle_quotient               ! Time -> angle

            !
            ! Precompute cosine and sine of the angles
            !

            cang = dcos(theta)
            sang = dsin(theta)

            !
            ! Construct transformation matrices
            !

            t_ir(1,1) = cang; t_ir(1,2) = -sang; t_ir(1,3) = 0.d0;
            t_ir(2,1) = sang; t_ir(2,2) = cang ; t_ir(2,3) = 0.d0;
            t_ir(3,1) = 0.d0; t_ir(3,2) = 0.d0 ; t_ir(3,3) = 1.d0;

            t_ir_dot(1,1) = -sang; t_ir_dot(1,2) = -cang; t_ir_dot(1,3) = 0.d0;
            t_ir_dot(2,1) = cang ; t_ir_dot(2,2) = -sang; t_ir_dot(2,3) = 0.d0;
            t_ir_dot(3,1:3) = 0.d0;

            state_out(1:3) = matmul(t_ir, state_in(1:3))
            state_out(4:6) = matmul(t_ir, state_in(4:6))

        end subroutine rotator

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Rotates a state in the inertial, ECLIPJ2000 (km, km/s) frame into the 
        !! CR3BP synodic frame (dimensionless, dimensionless)
        !
        ! REVISION HISTORY:
        ! 21 Mar 2018 - Initial version
        !
        !> @param[in] state_in, t
        !> @param[out] state_out
        !--------------------------------------------------------------------------- 

        subroutine synodic_rotate(state_in, t, state_out)

            use constants

            double precision, intent(in)   :: state_in(6)                                       ! Input state
            double precision, intent(in)   :: t                                                 ! Input time (ephemeris seconds)

            double precision, intent(out)  :: state_out(6)                                      ! Output state

            double precision               :: r(3), rdot(3)                                     ! Position, velocity
            double precision, parameter    :: theta_0 = 100.3762 * 3.141592 / 180.d0            ! Angle of the Earth at J2000

            double precision               :: t_ir(3, 3)                                        ! Transformation matrix
            double precision               :: t_irdot(3, 3)                                     ! Derivative of the transformation matrix

            double precision               :: total_angle                                       ! Total angle to rotate (inertial -> synodic)
            double precision               :: cang                                              ! Cosine of total_angle
            double precision               :: sang                                              ! Sine of total_angle
            
            !
            ! Get total angle from time to rotate
            !

            total_angle = theta_0 + (t * time_to_angle_quotient)

            !
            ! Pre-compute cos/sin of angle
            !

            cang = cos(total_angle)
            sang = sin(total_angle)

            !
            ! Separate state into position, velocity
            !

            r = state_in(1:3) * position_non_dimensionalise_quotient
            rdot = state_in(4:6) * velocity_non_dimensionalise_quotient

            !
            ! Construct transformation matrix
            !

            t_ir(1,1) = cang; t_ir(1, 2) = sang; t_ir(1,3) = 0.0
            t_ir(2,1) = -sang; t_ir(2,2) = cang; t_ir(2,3) = 0.0
            t_ir(3,1) = 0.0; t_ir(3,2) = 0.0; t_ir(3,3) = 1.0

            !
            ! Construct derivative of transformation matrix
            !

            t_irdot(1,1) = -sang; t_irdot(1,2) = cang; t_irdot(1,3) = 0.0
            t_irdot(2,1) = -cang; t_irdot(2,2) = -sang; t_irdot(2,3) = 0.0
            t_irdot(3,1) = 0.0; t_irdot(3,2) = 0.0; t_irdot(3,3) = 0.0

            !
            ! Multiply
            !

            state_out(1:3) = matmul(t_ir, r)
            state_out(4:6) = matmul(t_ir, rdot) + matmul(t_irdot, r)

            !
            ! Shift barycenter
            !

            state_out(1) = state_out(1) - mu3bp

        end subroutine synodic_rotate

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Rotates a state from the CR3BP synodic frame (dimensionless, dimensionless)
        !! into the inertial J2000 frame (km, km/s)
        !
        ! REVISION HISTORY:
        ! 21 Mar 2018 - Initial version
        !
        !> @param[in] state_in, t
        !> @param[out] state_out
        !--------------------------------------------------------------------------- 

        subroutine global_rotate(state_in, t, state_out)

            use constants

            double precision, intent(in)   :: state_in(6)                                       ! Input state
            double precision, intent(in)   :: t                                                 ! Input time (ephemeris seconds)

            double precision, intent(out)  :: state_out(6)                                      ! Output state

            double precision               :: r(3), rdot(3)                                     ! Position, velocity
            double precision, parameter    :: theta_0 = 100.3762 * 3.141592 / 180.d0            ! Angle of the Earth at J2000

            double precision               :: t_ir(3, 3)                                        ! Transformation matrix
            double precision               :: t_irdot(3, 3)                                     ! Derivative of the transformation matrix

            double precision               :: total_angle                                       ! Total angle to rotate (inertial -> synodic)
            double precision               :: cang                                              ! Cosine of total_angle
            double precision               :: sang                                              ! Sine of total_angle

            double precision               :: state_translate(6)                                ! State with translated barycentre
            
            ! Shift barycenter

            state_translate = state_in
            state_translate(1) = state_in(1) + mu3bp

            !
            ! Get total angle from time to rotate
            !

            total_angle = theta_0 + (t * time_to_angle_quotient)

            !
            ! Pre-compute cos/sin of angle
            !

            cang = cos(total_angle)
            sang = sin(total_angle)

            !
            ! Separate state into position, velocity
            !

            r = state_translate(1:3)
            rdot = state_translate(4:6)

            !
            ! Construct transformation matrix
            !

            t_ir(1,1) = cang; t_ir(1, 2) = sang; t_ir(1,3) = 0.0
            t_ir(2,1) = -sang; t_ir(2,2) = cang; t_ir(2,3) = 0.0
            t_ir(3,1) = 0.0; t_ir(3,2) = 0.0; t_ir(3,3) = 1.0

            !
            ! Construct derivative of transformation matrix
            !

            t_irdot(1,1) = -sang; t_irdot(1,2) = cang; t_irdot(1,3) = 0.0
            t_irdot(2,1) = -cang; t_irdot(2,2) = -sang; t_irdot(2,3) = 0.0
            t_irdot(3,1) = 0.0; t_irdot(3,2) = 0.0; t_irdot(3,3) = 0.0

            !
            ! Transpose and multiply
            !

            t_ir = transpose(t_ir)
            t_irdot = transpose(t_irdot)

            state_out(1:3) = matmul(t_ir, r)
            state_out(4:6) = matmul(t_ir, rdot) + matmul(t_irdot, r)

            !
            ! Dimensionalise
            !

            state_out(1:3) = state_out(1:3) * position_dimensionalise_quotient
            state_out(4:6) = state_out(4:6) * velocity_dimensionalise_quotient

        end subroutine global_rotate

end module
