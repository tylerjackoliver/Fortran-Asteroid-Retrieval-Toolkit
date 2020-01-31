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

module state_determination

    implicit none

        contains

            subroutine get_state(targ_can, time_lower, time_upper)

                ! //////////////////////////////////////////////////////////
                !
                ! Applies the state determination methodology from 
                ! Sanchez et. al.
                !
                ! TODO: Modularization
                !
                ! /////////////////////////////////////////////////////////

                use precision_kinds
                use constants

                implicit none

                double precision               :: epoch                                    ! Epoch: start time
                double precision               :: epoch_counter                            ! Epoch loop variable
                double precision               :: ang_can                                  ! Angle of the candidate with respect to the Sun
                double precision               :: tau_can                                  ! Period of the candidate
                double precision               :: epoch_upper                              ! Upper bound for search space
                double precision               :: last_checked                             ! Time candidate was last in pi/8
                double precision               :: dum                                      ! Dummy variable
                double precision               :: max_distance                             ! Maximum distance between the candidate and the Earth
                double precision, dimension(6) :: state_can, state_syn          ! States of the can.,Earth, and the synodic state of the candidate respectively

                CHARACTER(len=6)            :: abcorr                                   ! Abberation correcton string: initialised to NONE
                CHARACTER(len=5)            :: obs                                      ! Observing body string
                CHARACTER(len=12)           :: coord                                    ! Co-ordinate system of reference
                CHARACTER(len=19)           :: epoch_str                                ! String of the lower epoch bound
                CHARACTER(len=19)           :: epoch_upper_str                          ! String of the upper epoch bound
                CHARACTER(len=7)            :: targ_ear                                 ! String of the target body (Earth)

                INTEGER                     :: istate = 0                               ! Success flag, iteration counter, max. array locator

                CHARACTER(*),  intent(in)   :: targ_can                                 ! String of the candidate SPK ID

                double precision, intent(out)  :: time_lower                               ! Lower bound of the solution search time
                double precision, intent(out)  :: time_upper                               ! Upper bound of the solution search time


                max_distance = 0                                                        ! Initialise to trivially low to force update
                
                ! Ephemeris options
                
                targ_ear        = 'Earth'
                abcorr          = 'NONE'
                obs             = 'Sun'
                coord           = 'ECLIPJ2000'
                epoch_str       = 'Jan 1, 2020 00:00'       
                epoch_upper_str = 'Jan 1, 2100 00:00'

                ! call FURNSH('../data/de414.bsp')                    
                ! call FURNSH('../data/naif0008.tls')
                ! call FURNSH('../data/'//targ_can//'.bsp')

                ! call welcomemessage(targ_can)

                ! Program execution: Load the necessary SPICE kernels
                ! Future optimisation could move these into a meta-kernel

                ! Get the lower bound of the epoch in terms of ephemeris seconds (et)

                call STR2ET(epoch_str, epoch)

                ! As above, but with the upper bound

                call STR2ET(epoch_upper_str, epoch_upper)

                ! Initialise 'trip' counter to be the start of the domain

                last_checked = epoch                                                    ! To prevent errors on first run

                ! Initialise the period to which we enforce outside of pi/8 (seconds)

                tau_can = 365.25 * 86400.d0                                             ! One synodic periodic, seconds

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

                ! Refine coarse solution using a finer grid & get time at which
                ! Earth and candidate are farthest away

                time_lower = last_checked
                time_upper = last_checked+tau_can

                ! call UNLOAD('../data/de414.bsp')                    
                ! call UNLOAD('../data/naif0008.tls')
                ! call UNLOAD('../data/'//targ_can//'.bsp')

            end subroutine get_state

            subroutine time_str(str_in, time_out)

                use precision_kinds

                CHARACTER(*), intent(in) :: str_in

                double precision, intent(out) :: time_out
                
                call STR2ET(str_in, time_out)

            end subroutine time_str

            subroutine CANDIDATE_POSITION(t_epoch, state_can)

                ! Propagates forward a state in purely two-body motion 
                ! using the SPICE/NAIF routine PROP2B, which implements
                ! a Universal Variables propagation technique
                ! This subroutine assumes that the time variables are
                ! inputted as ephemeris time

                use precision_kinds
                use problem_parameters ! Gives access to targ_can

                double precision, intent(in)   :: t_epoch                               ! Initial time
                double precision, intent(out)  :: state_can(6)                          ! Final state

                CHARACTER(len=6)            :: abcorr                                   ! Abberation correcton string: initialised to NONE
                CHARACTER(len=5)            :: obs                                      ! Observing body string
                CHARACTER(len=12)           :: coord    
                double precision               :: dum

                ! call FURNSH('../data/de414.bsp')                    
                ! call FURNSH('../data/naif0008.tls')
                ! call FURNSH('../data/'//targ_can//'.bsp')

                abcorr = 'NONE'
                obs    = 'Sun'
                coord  = 'ECLIPJ2000'

                ! Get position from ephemeris

                call SPKEZR(targ_can, t_epoch, 'ECLIPJ2000', 'NONE', &
                            'Sun', state_can, dum)

                ! call UNLOAD('../data/de414.bsp')                    
                ! call UNLOAD('../data/naif0008.tls')
                ! call UNLOAD('../data/'//targ_can//'.bsp')
                
            end subroutine CANDIDATE_POSITION

            subroutine ROTATOR(state_in, epoch, tt, state_out)

                use precision_kinds
                use constants

                double precision, intent(in)   :: state_in(6)
                double precision, intent(in)   :: epoch
                double precision, intent(in)   :: tt
                double precision, intent(out)  :: state_out(6) 

                double precision               :: cang
                double precision               :: sang

                double precision               :: theta
                double precision               :: theta_dot
                double precision               :: t_ir(3, 3)
                double precision               :: t_ir_dot(3, 3)

                theta_dot = 2.d0 * pi / (86400.d0 * 365.25d0)
                theta     = (tt + epoch) * theta_dot

                cang = dcos(theta)
                sang = dsin(theta)

                t_ir(1,1) = cang; t_ir(1,2) = -sang; t_ir(1,3) = 0.d0;
                t_ir(2,1) = sang; t_ir(2,2) = cang ; t_ir(2,3) = 0.d0;
                t_ir(3,1) = 0.d0; t_ir(3,2) = 0.d0 ; t_ir(3,3) = 1.d0;

                t_ir_dot(1,1) = -sang; t_ir_dot(1,2) = -cang; t_ir_dot(1,3) = 0.d0;
                t_ir_dot(2,1) = cang ; t_ir_dot(2,2) = -sang; t_ir_dot(2,3) = 0.d0;
                t_ir_dot(3,1:3) = 0.d0;

                state_out(1:3) = matmul(t_ir, state_in(1:3))
                state_out(4:6) = matmul(t_ir, state_in(4:6))

            end subroutine ROTATOR

            subroutine SYNODIC_ROTATE(state_in, t, state_out)

                ! //////////////////////////////////////////////////////
                !
                ! Rotates state from the inertial frame (km, km/s) into
                ! the synodic frame (dimensionless)
                !
                ! //////////////////////////////////////////////////////

                use precision_kinds
                use constants

                double precision, intent(in)   :: state_in(6)
                double precision, intent(in)   :: t

                double precision, intent(out)  :: state_out(6)

                double precision               :: r(3), rdot(3)
                double precision, parameter    :: theta_0 = 100.3762 * 3.141592 / 180.d0           ! Angle of the Earth at J2000

                double precision               :: t_ir(3, 3)
                double precision               :: t_irdot(3, 3)

                double precision               :: total_angle
                double precision               :: cang
                double precision               :: sang 
                
                double precision, parameter    :: mu3bp = 3.003458d-06

                ! Get total angle - OK to make large as cos/sin functions will modulo 2pi the answer

                total_angle = theta_0 + (t * 2.d0 * pi) / (86400.d0 * 365.25d0)

                ! Pre-compute cos/sin of angle

                cang = cos(total_angle)
                sang = sin(total_angle)

                ! Non-dimensionalise first

                r = state_in(1:3) / au
                rdot = state_in(4:6) / au * 2 * pi / 365.25

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

    end module
