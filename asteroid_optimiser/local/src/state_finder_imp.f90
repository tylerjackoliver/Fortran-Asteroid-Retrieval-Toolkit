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
!   Areas for improvement
!   ~~~~~~~~~~~~~~~~~~~~~
!   - Define a meta-kernel to use to load invariant
!     kernels
!   - Computing the same thing (candidate/Ear state)
!     twice: is the memory tradeoff worth it?
!   - The code uses the last guess for the candidate period;
!     this is variant and only an approximation
!   - Parameterise by input body
!

module state_determination

    implicit none

        contains

            subroutine state_finder(targ_can,state_out,out_time,time_lower,time_upper)

                use precision_kinds                                         ! Defines double and quadruple prec. better
                use constants                                               ! Defines useful constants (pi, au, gm etc.)

                implicit none

                ! Variable declaration

                CHARACTER(len=*),intent(in) :: targ_can
                
                REAL(kind=dp),   intent(out) :: state_out(6)
                REAL(kind=dp),   intent(out) :: out_time
                REAL(kind=dp),   intent(out) :: time_lower
                REAL(kind=dp),   intent(out) :: time_upper

                REAL(kind=dp)                :: start 
                REAL(kind=dp)                :: finish

                call CPU_TIME(start)

                ! Load invariant kernels for SPICE

                call FURNSH('../data/de414.bsp')                    

                call FURNSH('../data/naif0008.tls')

                call get_state(targ_can,state_out,out_time,time_lower,time_upper)

                call CPU_TIME(finish)

                ! write(*,*) 'Done. Time to completion: ', finish-start, ' seconds.'
                ! write(*,*) "State_out is", state_out

                call UNLOAD('../data/de414.bsp')
                call UNLOAD('../data/naif0008.tls')
    
            end subroutine state_finder


            ! subroutine welcomemessage(targ_can)

            !     CHARACTER(*), intent(in) :: targ_can

            !     write(*,*) "    ______           __                      ___         __                  _     __   ____       __       _                 __   ______            __    ___________    ____  _______"
            !     write(*,*) "   / ____/___  _____/ /__________ _____     /   |  _____/ /____  _________  (_)___/ /  / __ \___  / /______(_)__ _   ______ _/ /  /_  __/___  ____  / /  _/_/ ____/   |  / __ \/_  __/ |"
            !     write(*,*) "  / /_  / __ \/ ___/ __/ ___/ __ `/ __ \   / /| | / ___/ __/ _ \/ ___/ __ \/ / __  /  / /_/ / _ \/ __/ ___/ / _ \ | / / __ `/ /    / / / __ \/ __ \/ /  / // /_  / /| | / /_/ / / /  / /"
            !     write(*,*) " / __/ / /_/ / /  / /_/ /  / /_/ / / / /  / ___ |(__  ) /_/  __/ /  / /_/ / / /_/ /  / _, _/  __/ /_/ /  / /  __/ |/ / /_/ / /    / / / /_/ / /_/ / /  / // __/ / ___ |/ _, _/ / /  / / "
            !     write(*,*) "/_/    \____/_/   \__/_/   \__,_/_/ /_/  /_/  |_/____/\__/\___/_/   \____/_/\__,_/  /_/ |_|\___/\__/_/  /_/\___/|___/\__,_/_/    /_/  \____/\____/_/  / //_/   /_/  |_/_/ |_| /_/ _/_/  "
            !     write(*,*) "                                                                                                                                                      |_|                        /_/    "

            ! end subroutine welcomemessage


            subroutine midmessage

                write(*,*) "Success! Refining..."

            end subroutine midmessage


            subroutine get_state(targ_can,state_out, time_count,time_lower,time_upper)

                use precision_kinds
                use constants

                REAL(kind=dp)               :: epoch                                    ! Epoch: start time
                REAL(kind=dp)               :: epoch_counter                            ! Epoch loop variable
                REAL(kind=dp)               :: ang_can                                  ! Angle of the candidate with respect to the Sun
                REAL(kind=dp)               :: ang_ear                                  ! Angle of the Earth with respect to the Sun
                REAL(kind=dp)               :: tau_can                                  ! Period of the candidate
                REAL(kind=dp)               :: a_can                                    ! Period/SMA of candidate
                REAL(kind=dp)               :: epoch_upper                              ! Upper bound for search space
                REAL(kind=dp)               :: last_checked                             ! Time candidate was last in pi/8
                REAL(kind=dp)               :: dum                                      ! Dummy variable
                REAL(kind=dp)               :: distance                                 ! Distance between the candidate and the Earth
                REAL(kind=dp)               :: max_distance                             ! Maximum distance between the candidate and the Earth
                REAL(kind=dp), dimension(6) :: state_can, state_ear                     ! States of the can. and Earth, respectively
                REAL(kind=dp), dimension(8) :: oe_can                                   ! OE vector

                CHARACTER(len=6)            :: abcorr                                   ! Abberation correcton string: initialised to NONE
                CHARACTER(len=5)            :: obs                                      ! Observing body string
                CHARACTER(len=12)           :: coord                                    ! Co-ordinate system of reference
                CHARACTER(len=19)           :: epoch_str                                ! String of the lower epoch bound
                CHARACTER(len=19)           :: epoch_upper_str                          ! String of the upper epoch bound
                CHARACTER(len=7)            :: targ_ear                                 ! String of the target body (Earth)

                INTEGER                     :: istate = 0                               ! Success flag, iteration counter, max. array locator

                CHARACTER(*),  intent(in)   :: targ_can                                 ! String of the candidate SPK ID

                REAL(kind=dp), intent(out)  :: state_out(6)                             ! State of the target at the reference epoch
                REAL(kind=dp), intent(out)  :: time_count                               ! Reference epoch for the target
                REAL(kind=dp), intent(out)  :: time_lower                               ! Lower bound of the solution search time
                REAL(kind=dp), intent(out)  :: time_upper                               ! Upper bound of the solution search time


                max_distance = 0                                                        ! Initialise to trivially low to force update
                targ_ear        = 'Earth'
                abcorr          = 'NONE'
                obs             = 'Sun'
                coord           = 'ECLIPJ2000'
                epoch_str       = 'Jan 1, 2020 00:00'       
                epoch_upper_str = 'Jan 1, 2050 00:00'

                ! call welcomemessage(targ_can)

                ! Program execution: Load the necessary SPICE kernels
                ! Future optimisation could move these into a meta-kernel

                call FURNSH(targ_can//'.bsp')                               ! Candidate ephemeris

                ! Get the epoch in terms of ephemeris seconds (et)

                call STR2ET(epoch_str, epoch)

                ! Get the upper bound for epoch

                call STR2ET(epoch_upper_str, epoch_upper)

                last_checked = epoch                                        ! To prevent errors on first run

                do epoch_counter = epoch, epoch_upper, 86400

                    ! Get the initial state of Earth and candidate

                    call SPKEZR(targ_can, epoch_counter, 'ECLIPJ2000', abcorr,&
                                obs, state_can, dum)
                    call SPKEZR(targ_ear, epoch_counter, 'ECLIPJ2000', abcorr,&
                                obs, state_ear, dum)

                    ! Get the orbital elements of the candidate (to compute orbital
                    ! period)

                    call OSCELT(state_can, epoch_counter, mu, oe_can)

                    ! Get the initial angle between Earth and candidate
                    ! Transcendental functions - very slow but no way around it

                    ang_can = datan2(state_can(2),state_can(1))
                    ang_ear = datan2(state_ear(2),state_ear(1))

                    a_can   = oe_can(1)/(1.0_dp - oe_can(2))                ! Calculate SMA from
                                                                            ! RP and e
                    tau_can = 2*pi*sqrt(a_can**3/mu)                        ! Period of the
                                                                            ! candidate, seconds

                    if (abs(ang_ear-ang_can) < pi/8) then

                        last_checked = epoch_counter

                    end if

                    if ((epoch_counter-last_checked) > tau_can) then

                        istate=1                                            ! Switch flag
                        exit                                                ! Exit the loop

                    end if

                end do

                ! Check if we found a solution

                if (istate /= 1) then

                    error stop "Solution not found. Stopping..."

                end if

                ! Now, let's continue...
                ! Variable values will still be saved, so let's just call them
                ! here

                ! call midmessage

                ! Allocate the size of the arrays based on the candidate

                !write(*,*) "Distance is: ", distance

                time_lower = last_checked
                time_upper = last_checked+tau_can

                do epoch_counter = time_lower, time_upper, 16600

                    ! Get the initial state of Earth and candidate

                    call SPKEZR(targ_can, epoch_counter, 'ECLIPJ2000', abcorr, obs,&
                                state_can,dum)

                    call SPKEZR(targ_ear, epoch_counter, 'ECLIPJ2000', abcorr, obs,&
                                state_ear,dum)

                    distance = sqrt((state_can(1)-state_ear(1))**2    &
                                            +(state_can(2)-state_ear(2))**2 &
                                            +(state_can(3)-state_ear(3))**2)
                    
                    if (distance .gt. max_distance) then

                        max_distance = distance
                        state_out = state_can
                        time_count = epoch_counter

                    end if

                end do

                ! Unload the SPICE kernel

                call UNLOAD(targ_can//'.bsp')

            end subroutine get_state

            subroutine prop_to_time(state_in, start_time, finish_time, state_out)

                ! Propagates forward a state in purely two-body motion 
                ! using the SPICE/NAIF routine PROP2B, which implements
                ! a Universal Variables propagation technique
                ! 
                ! The propagation time (finish_time-start_time) is split
                ! into periods of half-days and incremented
                !
                ! This subroutine assumes that the time variables are
                ! inputted as ephemeris time

                use precision_kinds
                use constants

                REAL(kind=dp), intent(in)   :: state_in(6)
                REAL(kind=dp), intent(in)   :: start_time
                REAL(kind=dp), intent(in)   :: finish_time

                REAL(kind=dp), intent(out)  :: state_out(6)

                REAL(kind=dp)               :: prop_time
                REAL(kind=dp)               :: prop_step = 33200.                 ! Propagation time-step: currently .5 day
                REAL(kind=dp)               :: prev_state(6)
                REAL(kind=dp)               :: remainder

                INTEGER                     :: num_of_iters
                INTEGER                     :: i

                ! Calculate propagation time (seconds)

                prop_time = finish_time - start_time    
                
                ! Calculate number of iterations

                num_of_iters = floor(prop_time/16600.)

                prev_state = state_in

                do i = 1,num_of_iters

                    call PROP2B(mu, prev_state, prop_step, state_out)
                    prev_state = state_out

                end do

                remainder = prop_time - num_of_iters * 16600.

                ! Propagate for the part 'missed' by the do-loop

                call PROP2B(mu, prev_state, remainder, state_out)

            end subroutine prop_to_time

            subroutine time_str(str_in, time_out)

                use precision_kinds

                CHARACTER(*), intent(in) :: str_in

                REAL(kind=dp), intent(out) :: time_out
                
                call STR2ET(str_in, time_out)

            end subroutine time_str

            subroutine CANDIDATE_POSITION(t_epoch,state_epoch,orig_state,state_can)

                use precision_kinds
                use constants
        
                REAL(kind=dp), intent(in)  :: t_epoch                               ! Initial time
                REAL(kind=dp), intent(in)  :: state_epoch                           ! Epoch for the desired state
                REAL(kind=dp), intent(in)  :: orig_state(6)                         ! Original state
                
                REAL(kind=dp), intent(out) :: state_can(6)                          ! Final state
        
                ! Call PROP2B with the propagation time being the difference
                ! between t_epoch and state_epoch
        
                call PROP2B(mu, orig_state, t_epoch - state_epoch, state_can)
        
            end subroutine CANDIDATE_POSITION

            subroutine ROTATOR(state_in, epoch, tt, state_out)

                ! ////////////////////////////////////////////////////////////////////////////
                !
                ! Applies relations from Sanchez et. al. to "rotate" a target state by means
                ! of altering its longitude
                !
                ! ///////////////////////////////////////////////////////////////////////////
        
                use precision_kinds
                use constants
        
                REAL(kind=dp), intent(in)   :: state_in(6)                                                                              ! Input state; pre-rotation
                REAL(kind=dp), intent(in)   :: epoch                                                                                    ! Epoch of the state
                REAL(kind=dp), intent(in)   :: tt                                                                                       ! Transfer time
                
                REAL(kind=dp), intent(out)  :: state_out(6)                                                                             ! Output state; post-rotation
        
                REAL(kind=dp)               :: elts(8)                                                                                  ! Elements array
        
                ! Convert from state -> orbital elements
                ! Epoch does nothing here! So still need to rotate
                ! by (epoch + tt) later on
        
                call OSCELT(state_in, epoch, mu, elts)
        
                ! Apply relations from Sanchez et. al. The following code
                ! block assumes that the time t is inputted as ephemeris
                ! seconds past J2000 (and thus that it is linked innately to
                ! the epoch of the states used in this work.
        
                if (abs(elts(4)) < 1.0d-6) then
        
                    elts(4) = 0
                    elts(5) = elts(5) + (epoch+tt)*(2.d0*pi)/(365.25d0*86400.d0)                                                                   ! Add on phasing equal to the transfer time                                           
        
                else
        
                    elts(4) = elts(4) + (epoch+tt)*(2.d0 * pi)/(365.25d0*86400.d0)                                                                   ! Add on phasing equal to the transfer time
                    elts(5) = 0
        
                end if
        
                ! Convert from orbital elements -> state at the later epoch
        
                call CONICS(elts, epoch + tt, state_out)
        
            end subroutine ROTATOR

    end module
