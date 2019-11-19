        program main

            use precision_kinds
            use constants
            
            implicit none

            REAL(kind=dp)           :: state(6),stateout(6)
            REAL(kind=dp)           :: t_ir(3,3), t_irdot(3,3)
            REAL(kind=dp)           :: rot_ang, cang, sang
            REAL(kind=sp)           :: dum

            INTEGER(kind=sp)        :: io_state
            INTEGER(kind=sp)        :: i = 0
            INTEGER(kind=sp)        :: unitno = 50
            INTEGER(kind=sp)        :: outunit= 51

            CHARACTER(len=21)       :: infile = 'output_data_mod.dat'
            CHARACTER(len=21)       :: output = 'output_data_fin.dat'

            ! Open file

            open(unit=unitno, file=infile)
            open(unit=outunit,file=output)
            ! Define rotation matrices

            rot_ang = 100.3762*pi/180.d0

            cang = dcos(rot_ang)
            sang = dsin(rot_ang)

            t_ir(1,1) = cang
            t_ir(2,1) = sang
            t_ir(3,1) = 0

            t_ir(1,2) = -sang
            t_ir(2,2) = cang
            t_ir(3,2) = 0

            t_ir(1,3) = 0
            t_ir(2,3) = 0
            t_ir(3,3) = 1


            t_irdot(1,1) = -sang
            t_irdot(2,1) = cang
            t_irdot(3,1) = 0

            t_irdot(1,2) = -cang
            t_irdot(2,2) = -sang
            t_irdot(3,2) = 0

            t_irdot(1:3,3) = 0


            do

                i = i+1

                read(unitno,*,iostat=io_state) dum, state(1), state(2),&
                    state(3), state(4), state(5), state(6)

                if (io_state /= 0) then

                    exit

                end if

                stateout(1) = stateout(1) + 3.0334d-06

                stateout(1:3) = matmul(t_ir,state(1:3))
                stateout(4:6) = matmul(t_ir,state(4:6)) + &
                    matmul(t_irdot,state(1:3))

                ! Dimensionalise the system

                stateout(1:3) = stateout(1:3)*au
                stateout(4:6) = &
                    stateout(4:6)*2.d0*pi*au/(86400.d0*365.25d0)

                write(outunit,'(3F16.3,3F16.8)') stateout(1:6)

                write(*,*) "Finished processing item", i

            end do

            close(unitno)
            close(outunit)

        end program main
