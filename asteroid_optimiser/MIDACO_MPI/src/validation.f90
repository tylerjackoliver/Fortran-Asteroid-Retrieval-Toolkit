module validation

    use precision_kinds
    use constants
    use tbp_ode

    implicit none

    double precision    :: integration_tol = 1.d-09

    contains

        subroutine least_squares_target(f, x, initial_state, desired_state)

            use integrators

            implicit none

            double precision, intent(out)   :: f(1)

            double precision, intent(in)    :: x(4)
            double precision, intent(in)    :: initial_state(6)
            double precision, intent(in)    :: desired_state(6)

            type(t_integrator)              :: integrator
            double precision                :: x0(6)
            double precision                :: xF(6)
            double precision                :: time 
        
            time = x(4)

            ! Initialise initial state

            x0 = (/initial_state(1), initial_state(2), initial_state(3), &
                   x(1), x(2), x(3)/)

            ! Integrate forward

            call PROPAGATE_STATE(x0, time, xF)

            ! Compute least-squares distance

            f = norm2(xF - desired_state)

        end subroutine least_squares_target

        subroutine validate_transfer(t0, tt, tmani, r1, v1, r2, v2, xtarg, J, validated_dv)

            use precision_kinds
            use constants

            implicit none

            double precision, intent(in)    :: t0
            double precision, intent(in)    :: tt
            double precision, intent(in)    :: tmani
            double precision, intent(in)    :: r1(6)
            double precision, intent(in)    :: v1(3)
            double precision, intent(in)    :: r2(6)
            double precision, intent(in)    :: v2(3)
            double precision, intent(in)    :: xtarg(6)
            
            integer,          intent(in)    :: J
 
            double precision, intent(out)   :: validated_dv

            double precision                :: MJD_r1
            double precision                :: MJD_r2
            double precision                :: MJD_xtarg

            double precision                :: lambert_begin(6)
            double precision                :: lambert_end(6)
            double precision                :: lambert_end_synodic(6)
            double precision                :: r1_non_dim(6)
            double precision                :: r1_synodic(6)
            double precision                :: r1_validated(6)
            double precision                :: r1_validated_non_dim(6)
            double precision                :: r1_validated_global(6)
            double precision                :: r2_non_dim(6)
            double precision                :: r2_global(6)

            double precision, external      :: unitim

            double precision                :: dv1
            double precision                :: dv2

            double precision                :: tt_synodic

            ! Load the CSPICE kernels

            call FURNSH('../data/naif0008.tls')

            ! Determine times

            MJD_r1 = UNITIM(t0, 'ET', 'JED') - 2400000.5                  ! t0
            MJD_r2 = UNITIM(t0+tt, 'ET', 'JED') - 2400000.d0              ! t0+tt
            MJD_xtarg = UNITIM(t0+tt+abs(tmani), 'ET', 'JED') ! Final time

            tt_synodic = tt * 2.d0 * pi / (365.25d0 * 86400.d0)

            ! Apply velocity v1 to the initial condition

            lambert_begin       = r1
            lambert_begin(4:6)  = r1(4:6) + v1

            ! Rotate into the synodic frame

            call SYNODIC_ROTATE(lambert_begin, MJD_r1, r1_synodic)

            ! Wiggle to get the correct state

            !print *, "I have r1_synodic as", r1_synodic, "and I'm about to give it a go"

            call GET_DESIRED_STATE(r1_synodic, r2, tt_synodic, r1_validated)

            ! Determine velocity requirement by dimensionalising
            
            call GLOBAL_ROTATE(r1_validated, MJD_r1, r1_validated_non_dim)

            r1_validated_global(1:3) = r1_validated_non_dim(1:3) * au
            r1_validated_global(4:6) = r1_validated_non_dim(4:6) * au * 2.d0 * pi / (365.25d0 * 86400.d0)

            dv1 = norm2(r1_validated_global(4:6) - r1(4:6))

            ! Propagate forward

            call PROPAGATE_STATE(r1_validated, tt_synodic, lambert_end_synodic)

            ! Move back into the global frame

            call GLOBAL_ROTATE(lambert_end_synodic, MJD_r2, lambert_end)

            ! Dimensionalise

            lambert_end(1:3) = lambert_end(1:3) * au
            lambert_end(4:6) = lambert_end(4:6) * au * 2.d0 * pi / (86400.d0 * 365.25d0)

            ! r2 is now synodic, so also rotate that to work out velocity requirements

            call GLOBAL_ROTATE(r2, MJD_r2, r2_non_dim)

            ! Dimensionalise

            r2_global(1:3) = r2_non_dim(1:3) * au
            r2_global(4:6) = r2_non_dim(4:6) * au * 2.d0 * pi / (365.25d0 * 86400.d0)

            ! Compute the new velocity

            validated_dv = dv1

        end subroutine validate_transfer

        subroutine PROPAGATE_STATE(initial, time, final)

            use integrators

            implicit none

            double precision, intent(in)    :: initial(6)
            double precision, intent(in)    :: time
            double precision, intent(out)    :: final(6)

            type(t_integrator)              :: integrator

            ! Set initial things for the integrator

            print *, "Tolerance: ", integration_tol

            call integrator%default_options()
            call integrator%set(tbp_eom, 6, 0)
            print *, "Time: ", time
            call integrator%options(RelTol=[integration_tol,integration_tol,&
            integration_tol,integration_tol,integration_tol,integration_tol],&
            AbsTol=[integration_tol,integration_tol,integration_tol,integration_tol,&
            integration_tol,integration_tol], verbose=.true.,&
            progress=.true.)    
            call integrator%initial_conditions(0.d0, initial)
            call integrator%integrate(time)
            final = integrator%y
            print *, final
            
        end subroutine PROPAGATE_STATE

        subroutine GET_DESIRED_STATE(initial, desired, time, initial_validated)

            double precision, intent(in)    :: initial(6)
            double precision, intent(in)    :: desired(6)
            double precision, intent(in)    :: time
            double precision, intent(out)    :: initial_validated(6)

            integer, parameter              :: O = 1            ! Num objective funcs
            integer, parameter              :: N = 4            ! Num params                                            
            integer, parameter              :: NI = 0           ! Number of integer constraints
            integer, parameter              :: M = 0            ! Number of contrains
            integer, parameter              :: ME = 0           ! Number of equality constraints
            integer, parameter              :: LIW = 5000       ! Integer workspace
            integer, parameter              :: LRW = 20000      ! Real workspace
            integer, parameter              :: LPF = 1000       ! Number of pareto front
        
            double precision                :: XL(4), XU(4)     ! Upper and lower bounds
            double precision                :: XOPT(4)          ! Optimisation variables
            double precision                :: F(1)             ! Objectives
            double precision                :: G(4)             ! Constraint arrays; initialised but never used
            double precision                :: PARAM(13)        ! MIDACO parameters
            double precision                :: RW(LRW), PF(LPF) ! Workspace and pareto front
        
            integer                         :: optim_flag       ! Optimiser information flag
            integer                         :: optim_stop       ! Optimiser stopping variable
            integer                         :: IW(LIW)          ! Integer workspace

            double precision                :: v1(2)
            double precision                :: v2(2)
            double precision                :: v3(2)

            integer                         :: max_time         ! Maximum walltime for optimisations
            integer                         :: max_eval         ! Maximum function evaluations
            integer                         :: print_eval       ! How often to print
            integer                         :: save_to_file     ! Output verbosity

            character(len=*), parameter     :: key = 'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]'

            print_eval = 1000000
            max_eval = 10000e3
            max_time = 5.d0 * 60

            PARAM( 1) = 0.0D0       ! ACCURACY
            PARAM( 2) = 0.0D0       ! SEED
            PARAM( 3) = 0.0D0       ! FSTOP
            PARAM( 4) = 0.0D0       ! ALGOSTOP
            PARAM( 5) = 500.D3      ! EVALSTOP
            PARAM( 6) = 0.D0        ! FOCUS
            PARAM( 7) = 0.0D0       ! ANTS
            PARAM( 8) = 0.0D0       ! KERNEL
            PARAM( 9) = 0.0D0       ! ORACLE
            PARAM(10) = 0.D0        ! PARETOMAX
            PARAM(11) = 0.0D0       ! EPSILON  
            PARAM(12) = 0.d0        ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0       ! CHARACTER 

            v1 = (/initial(4) * .5, initial(4) * 1.5/)
            v2 = (/initial(5) * .5, initial(5) * 1.5/)
            v3 = (/initial(6) * .5, initial(6) * 1.5/)

            XL = (/minval(v1), minval(v2), minval(v3), time * .8/)
            XU = (/maxval(v1), maxval(v2), maxval(v3), time * 1.2/)

            xopt = xl

            F(1) = 1000.d0

            call midaco_print(1, print_eval, save_to_file, optim_flag, optim_stop, F, G, XOPT, &
            XL, XU, O, N, NI, M, ME, RW, PF, max_eval, max_time, param, 1, 0, key)

            do while (F(1) .gt. 1.d-10 .and. optim_stop .eq. 0) 

                ! Evaluate objective function and constraints (none)

                call least_squares_target(f, xopt, initial, desired)

                ! Call MIDACO

                call midaco(1, O, N, NI, M, ME, XOPT, F, G, XL, XU, optim_flag, optim_stop, &
                        param, rw, lrw, iw, liw, pf, lpf, key)

                ! Print again

                call midaco_print(2, print_eval, save_to_file, optim_flag, optim_stop, F, G, XOPT, &
                XL, XU, O, N, NI, M, ME, RW, PF, max_eval, max_time, param, 1, 0, key)

            end do

            initial_validated = (/initial(1), initial(2), initial(3), xopt(1), xopt(2), xopt(3)/)

        end subroutine GET_DESIRED_STATE

end module validation