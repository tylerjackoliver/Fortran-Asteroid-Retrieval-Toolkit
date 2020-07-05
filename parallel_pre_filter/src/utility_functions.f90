module utility_functions

    use constants
    use precision_kinds

    implicit none

    contains

        subroutine global_rotate(state_in, t, state_out)

            ! //////////////////////////////////////////////////////
            !
            ! Rotates state from the synodic frame (dimensionless, dimensionless)
            ! into the global frame (dimensionless, dimensionless)
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
            
            double precision, parameter    :: mu3bp = 3.0032080443d-06

            ! Get total angle - OK to make large as cos/sin functions will modulo 2pi the answer

            total_angle = theta_0 + (t * 2.d0 * pi) / (86400.d0 * 365.25d0)

            ! Pre-compute cos/sin of angle

            cang = cos(total_angle)
            sang = sin(total_angle)

            ! Non-dimensionalise first

            r = state_in(1:3)
            rdot = state_in(4:6)

            !
            ! Shift barycentre
            !

            r(1) = state_in(1) + mu3bp

            ! Construct T

            t_ir(1,1) = cang;   t_ir(1, 2) = sang;  t_ir(1,3) = 0.0
            t_ir(2,1) = -sang;  t_ir(2,2) = cang;   t_ir(2,3) = 0.0
            t_ir(3,1) = 0.0;    t_ir(3,2) = 0.0;    t_ir(3,3) = 1.0

            ! Construct T_dot

            t_irdot(1,1) = -sang;   t_irdot(1,2) = cang;    t_irdot(1,3) = 0.0
            t_irdot(2,1) = -cang;   t_irdot(2,2) = -sang;   t_irdot(2,3) = 0.0
            t_irdot(3,1) = 0.0;     t_irdot(3,2) = 0.0;     t_irdot(3,3) = 0.0

            !
            ! Transpose (since this is the inertial -> synodic matrix) and multiply
            !

            t_ir = transpose(t_ir)
            t_irdot = transpose(t_irdot)

            state_out(1:3) = matmul(t_ir, r)
            state_out(4:6) = matmul(t_ir, rdot) + matmul(t_irdot, r)

        end subroutine global_rotate

end module utility_functions
