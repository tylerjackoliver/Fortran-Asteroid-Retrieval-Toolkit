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
!   DynamicalArrays (included as a module)
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
!


program main

    use omp_lib
    implicit none

    integer :: nthreads = 8 ! Number of threads to use
    integer :: unitno, candidatecount, itercount, threadno
    integer :: i, j, k ! Loop variables
    integer :: bodyID

    real :: a_b, e_b, inb, rp_t, e_t, in_t, start, finish
    real :: mindV, au
    real :: t1, t2, t3, t4, t5, date, o, Om, M, p


    integer, dimension(17719) :: database

    real, dimension(1353, 3) :: bodyData
    real, dimension(3046544,3) :: targetData

    character(len=260) :: bodystring, filepath
    character(len=8) :: fmt = '(I7.7)'

    candidatecount = 0
    itercount = 0
    au =  1.495978707e8
    filepath = 'corrected_oe/'

    ! Print welcome message

    call system('clear')

    write(*,*) "============================================================"
    write(*,*) " "
    write(*,*) "ASTEROID PRUNING PRE-FILTER: FORTRAN08 IMPLEMENTATION"
    write(*,*) " "
    write(*,*) "============================================================"

    ! Initialise OMP and delete current results.txt

    call OMP_SET_NUM_THREADS(nthreads)    
    call file_check()

    ! Load data sets

    write(*,*) "Loading data sets..."

    open(unit=3, file="target_csv.csv")    ! Target OEs

    do i = 1,3046544

        read(3,*) targetData(i,1), targetData(i,2), targetData(i,3),t1,t2,t3,t4,t5

    end do

    close(3)

    open(unit=20, file="Database_SPKs.csv")   ! Database file
    
    read(20,*) ! Read once to avoid 'spk_id' header

    do i = 1,size(database)
	    read(20,*) database(i)
    end do

    close(20)

    ! Write progress message and begin timing

    write(*,*) "Loading complete. Beginning compute sequence."
    start = omp_get_wtime()

    !$OMP PARALLEL DO PRIVATE(I,J,K, bodyid, bodystring, unitno, mindV, bodyData, a_b, e_b, inb, rp_t, e_t, in_t)
    do i = 1, size(database)! All bodies

        bodyID = database(i)
        write(bodystring, fmt) bodyID

        unitno = OMP_GET_THREAD_NUM() + 30

        open(unit=unitno, file=trim(filepath)//trim(bodystring)//'_OE_corr.csv')
    
        do j = 1,1353
            read(unitno,*) date, bodyData(j,1), bodyData(j,2), bodyData(j,3), o, O, M, p
        end do

        close(unitno)

        exitloop: do j = 272,1353
        ! exitloop:  do j = 220,220


            a_b = bodyData(j,1)*au
            e_b = bodyData(j,2)
            inb = bodyData(j,3)

            do k = 1,3046544
                
                rp_t = targetData(k,1)*au
                e_t = targetData(k,2)
                in_t = targetData(k,3)

                call cheapest(a_b, e_b, inb, rp_t, e_t, in_t, mindV)

                if (mindV < .7) then 

                    call file_add(bodyID, unitno, mindv)
                    candidatecount = candidatecount + 1
		            exit exitloop

                end if

            end do

        end do exitloop

        itercount = itercount + 1

        call system('clear')
        write(*,*) "============================================================"
        write(*,*) " "
        write(*,*) "ASTEROID PRUNING PRE-FILTER: FORTRAN08 IMPLEMENTATION"
        write(*,*) " "
        write(*,*) "============================================================"
        write(*,*) "Completed processing object number ", itercount
        write(*,*) "Candidates found: ", candidatecount

    end do
    !$OMP END PARALLEL DO

    call system('clear')
    write(*,*) "============================================================"
    write(*,*) " "
    write(*,*) "ASTEROID PRUNING PRE-FILTER: FORTRAN08 IMPLEMENTATION"
    write(*,*) " "
    write(*,*) "============================================================"

    finish = omp_get_wtime()

    write(*,*) "==================== Scan completed ========================"
    write(*,*) "Objects scanned: ", size(database)
    write(*,*) "Candidates found: ", candidatecount
    write(*,*) "Compute time: ", finish-start, "seconds."

end program main


subroutine file_add(bodyID, unitnum, mindv)

    ! Subroutine intents

    character(len=13) :: filepath ! Define address of file to be read
    integer, intent(in) :: bodyID
    integer, intent(in) :: unitnum
    real, intent(in) :: mindv
    logical :: exist

    filepath = "results.txt"

    inquire(file=filepath, exist=exist)

    if (exist) then

      open(unitnum, file=filepath, status="old", position="append", action="write")

    else

      open(unitnum, file=filepath, status="new", action="write")

    end if

    write(unitnum, *) bodyID, mindv
    close(unitnum)

end subroutine file_add


subroutine file_check()

    character(len=13) :: filepath
    logical :: exist

    filepath = 'results.txt'

    inquire(file=filepath, exist=exist)

    if (exist) then

        open(20, file=filepath)
        close(20, status='delete')

    end if

end subroutine file_check


subroutine cheaper(a_b, e_b, inb, rp_t, e_t, in_t, cheap)

    real, intent(in) :: a_b, e_b, inb, rp_t, e_t, in_t
    real, intent(out) :: cheap
    real :: mu_s, a_t, ra_t, rp_b, ra_b, a_int, in_b
    real :: dV1, dV2, dV3, dV4, dV5, dV6, dV7, dV8
    real :: dVi1, dVi2, dVi3, dVi4, dVi5, dVi6, dVi7, dVi8, dVi9, dVi10, dVi11, dVi12, dVi13, dVi14, dVi15, dVi16
    real :: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16
    double precision :: pi = 4.*datan(1.d0)

    mu_s = 1.327124400189e11

    in_b = inb*pi/180

    ra_b = a_b*(1+e_b)
    rp_b = a_b*(1-e_b)

    a_t = rp_t/(1-e_t) 
    ra_t = a_t*(1+e_t) 

    ! There are 16 cases of different transfer possibilities. 

    !! Cases 1-4: target rp_t, perform periapsis changing apoapsis first

    a_int = .5*(rp_b+rp_t) 

    dV1 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    dV2 = sqrt(mu_s*(2/rp_t-1/a_t)) - sqrt(mu_s*(2/rp_t-1/a_int)) 

    dVi1 = 2*sqrt((mu_s/a_b)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 

    if (rp_t > rp_b) then

        dVi2 = 2*sqrt((mu_s/a_int) * (rp_t/rp_b))*sin(abs(in_b-in_t)/2) 
        dVi3 = 2*sqrt((mu_s/a_int) * (rp_b/rp_t))*sin(abs(in_b-in_t)/2) 
	
    else

        dVi2 = 2*sqrt((mu_s/a_int)*(rp_b/rp_t))*sin(abs(in_b-in_t)/2) 
        dVi3 = 2*sqrt((mu_s/a_int)*(rp_t/rp_b))*sin(abs(in_b-in_t)/2) 
	
    end if

    dVi4 = 2*sqrt((mu_s/a_t)*(rp_b/ra_t))*sin(abs(in_b-in_t)/2) 

    v1 = sqrt(dV1**2+dVi1**2) + sqrt(dV2**2) 
    v2 = sqrt(dV1**2+dVi2**2) + sqrt(dV2**2) 
    v3 = sqrt(dV1**2) + sqrt(dVi3**2+dV2**2) 
    v4 = sqrt(dV1**2) + sqrt(dVi4**2 + dV2**2) 

    !! Cases 5:8: target rp_t, perform apoapsis changing periapsis first

    a_int = .5*(rp_t+ra_b) 

    dV3 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b)) 
    dV4 = sqrt(mu_s*(2/rp_t-1/a_t))-sqrt(mu_s*(2/rp_t-1/a_int)) 

    dVi5 = 2*sqrt(mu_s/a_b*(rp_b/ra_b))*sin(abs(in_b-in_t)/2) 

    if (ra_b > rp_t) then

        dVi6 = 2*sqrt((mu_s/a_int)*(rp_t/ra_b))*sin(abs(in_b-in_t)/2) 
        dVi7 = 2*sqrt((mu_s/a_int)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 

    else
        
        dVi6 = 2*sqrt((mu_s/a_int)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 
        dVi7 = 2*sqrt((mu_s/a_int)*(rp_t/ra_b))*sin(abs(in_b-in_t)/2) 
        
    end if

    dVi8 = 2*sqrt((mu_s/a_t)*ra_t/rp_t)*sin(abs(in_b-in_t)/2) 

    v5 = sqrt(dVi5**2+dV3**2)+sqrt(dV4**2) 
    v6 = sqrt(dV3**2+dVi6**2)+sqrt(dV4**2) 
    v7 = sqrt(dV3**2) + sqrt(dV4**2+dVi7**2) 
    v8 = sqrt(dV3**2) + sqrt(dV4**2+dVi8**2) 

    !! Cases 9-12: target ra_t, perform periapsis changing apoapsis first

    a_int = .5*(ra_t+rp_b) 

    dV5 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    dV6 = sqrt(mu_s*(2/ra_t-1/a_t)) - sqrt(mu_s*(2/ra_t-1/a_int)) 

    dVi9 = 2*sqrt((mu_s/a_b)*(ra_b/rp_b))*sin(abs(in_b-in_t)/2) 

    if (ra_t > rp_b) then

        dVi10 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
        dVi11 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
        
    else

        dVi10 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
        dVi11 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
        
    end if

    dVi12 = 2*sqrt((mu_s/a_t)*(rp_t)/ra_t)*sin(abs(in_b-in_t)/2) 

    v9 = sqrt(dVi9**2+dV5**2)+sqrt(dV6**2) 
    v10 = sqrt(dV5**2+dVi10**2)+sqrt(dV6**2) 
    v11 = sqrt(dV5**2) + sqrt(dV6**2+dVi11**2) 
    v12 = sqrt(dV5**2) + sqrt(dV6**2+dVi12**2) 

    !! Cases 13-16: target ra_t, perform ap changing peri first

    a_int = .5*(ra_b+ra_t) 

    dV7 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b)) 
    dV8 = sqrt(mu_s*(2/ra_t-1/a_t))-sqrt(mu_s*(2/ra_t-1/a_int)) 

    dVi13 = 2*sqrt(mu_s/a_b*(rp_b/ra_b))*sin(abs(in_b-in_t)) 

    if (ra_b > ra_t) then

        dVi14 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 
        dVi15 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 

    else

        dVi14 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 
        dVi15 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 

    end if

    dVi16 = 2*sqrt(mu_s/a_int*(ra_t/rp_t))*sin(abs(in_b-in_t)) 

    v13 = sqrt(dVi13**2+dV7**2)+sqrt(dV8**2) 
    v14 = sqrt(dV7**2+dVi14**2)+sqrt(dV8**2) 
    v15 = sqrt(dV7**2) + sqrt(dV8**2+dVi15**2) 
    v16 = sqrt(dV7**2) + sqrt(dV8**2+dVi16**2) 

    cheap = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16)

end subroutine cheaper

subroutine cheapest(a_b, e_b, inb, rp_t, e_t, in_t, cheap)

    !
    ! Computes the minimum transfer velocity between a target 
    !

    real, intent(in)    :: a_b, e_b, inb, rp_t, e_t, in_t
    real, intent(out)   :: cheap
    
    real                :: mu_s, a_t, ra_t, rp_b, ra_b, a_int, in_b
    real                :: dV1, dV2, dV3, dV4, dV5, dV6, dV7, dV8
    real                :: dVi1, dVi2, dVi3, dVi4, dVi5, dVi6, dVi7, dVi8, dVi9, dVi10, dVi11, dVi12, dVi13, dVi14, dVi15, dVi16
    real                :: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16
    
    double precision    :: pi = 4.*datan(1.d0)

    mu_s = 1.327124400189e11

    in_b = inb*pi/180

    ra_b = a_b*(1+e_b)
    rp_b = a_b*(1-e_b)

    a_t = rp_t/(1-e_t) 
    ra_t = a_t*(1+e_t) 

    ! There are 16 cases of different transfer possibilities. 

    !! Cases 1-4: target rp_t, perform periapsis changing apoapsis first

    a_int = .5*(rp_b+rp_t) 

    dV1 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    dV2 = sqrt(mu_s*(2/rp_t-1/a_t)) - sqrt(mu_s*(2/rp_t-1/a_int)) 

    dVi1 = 2*sqrt((mu_s/a_b)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 

    if (rp_t > rp_b) then

        dVi2 = 2*sqrt((mu_s/a_int) * (rp_t/rp_b))*sin(abs(in_b-in_t)/2) 
        dVi3 = 2*sqrt((mu_s/a_int) * (rp_b/rp_t))*sin(abs(in_b-in_t)/2) 
	
    else

        dVi2 = 2*sqrt((mu_s/a_int)*(rp_b/rp_t))*sin(abs(in_b-in_t)/2) 
        dVi3 = 2*sqrt((mu_s/a_int)*(rp_t/rp_b))*sin(abs(in_b-in_t)/2) 
	
    end if

    dVi4 = 2*sqrt((mu_s/a_t)*(rp_b/ra_t))*sin(abs(in_b-in_t)/2) 

    v1 = sqrt(dV1**2+dVi1**2) + sqrt(dV2**2) 
    v2 = sqrt(dV1**2+dVi2**2) + sqrt(dV2**2) 
    v3 = sqrt(dV1**2) + sqrt(dVi3**2+dV2**2) 
    v4 = sqrt(dV1**2) + sqrt(dVi4**2 + dV2**2) 

    !! Cases 5:8: target rp_t, perform apoapsis changing periapsis first

    a_int = .5*(rp_t+ra_b) 

    dV3 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b)) 
    dV4 = sqrt(mu_s*(2/rp_t-1/a_t))-sqrt(mu_s*(2/rp_t-1/a_int)) 

    dVi5 = 2*sqrt(mu_s/a_b*(rp_b/ra_b))*sin(abs(in_b-in_t)/2) 

    if (ra_b > rp_t) then

        dVi6 = 2*sqrt((mu_s/a_int)*(rp_t/ra_b))*sin(abs(in_b-in_t)/2) 
        dVi7 = 2*sqrt((mu_s/a_int)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 

    else
        
        dVi6 = 2*sqrt((mu_s/a_int)*(ra_b/rp_t))*sin(abs(in_b-in_t)/2) 
        dVi7 = 2*sqrt((mu_s/a_int)*(rp_t/ra_b))*sin(abs(in_b-in_t)/2) 
        
    end if

    dVi8 = 2*sqrt((mu_s/a_t)*ra_t/rp_t)*sin(abs(in_b-in_t)/2) 

    v5 = sqrt(dVi5**2+dV3**2)+sqrt(dV4**2) 
    v6 = sqrt(dV3**2+dVi6**2)+sqrt(dV4**2) 
    v7 = sqrt(dV3**2) + sqrt(dV4**2+dVi7**2) 
    v8 = sqrt(dV3**2) + sqrt(dV4**2+dVi8**2) 

    !! Cases 9-12: target ra_t, perform periapsis changing apoapsis first

    a_int = .5*(ra_t+rp_b) 

    dV5 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    dV6 = sqrt(mu_s*(2/ra_t-1/a_t)) - sqrt(mu_s*(2/ra_t-1/a_int)) 

    dVi9 = 2*sqrt((mu_s/a_b)*(ra_b/rp_b))*sin(abs(in_b-in_t)/2) 

    if (ra_t > rp_b) then

        dVi10 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
        dVi11 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
        
    else

        dVi10 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
        dVi11 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
        
    end if

    dVi12 = 2*sqrt((mu_s/a_t)*(rp_t)/ra_t)*sin(abs(in_b-in_t)/2) 

    v9 = sqrt(dVi9**2+dV5**2)+sqrt(dV6**2) 
    v10 = sqrt(dV5**2+dVi10**2)+sqrt(dV6**2) 
    v11 = sqrt(dV5**2) + sqrt(dV6**2+dVi11**2) 
    v12 = sqrt(dV5**2) + sqrt(dV6**2+dVi12**2) 

    !! Cases 13-16: target ra_t, perform ap changing peri first

    a_int = .5*(ra_b+ra_t) 

    dV7 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b)) 
    dV8 = sqrt(mu_s*(2/ra_t-1/a_t))-sqrt(mu_s*(2/ra_t-1/a_int)) 

    dVi13 = 2*sqrt(mu_s/a_b*(rp_b/ra_b))*sin(abs(in_b-in_t)) 

    if (ra_b > ra_t) then

        dVi14 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 
        dVi15 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 

    else

        dVi14 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 
        dVi15 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 

    end if

    dVi16 = 2*sqrt(mu_s/a_int*(ra_t/rp_t))*sin(abs(in_b-in_t)) 

    v13 = sqrt(dVi13**2+dV7**2)+sqrt(dV8**2) 
    v14 = sqrt(dV7**2+dVi14**2)+sqrt(dV8**2) 
    v15 = sqrt(dV7**2) + sqrt(dV8**2+dVi15**2) 
    v16 = sqrt(dV7**2) + sqrt(dV8**2+dVi16**2) 

    cheap = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16)

end subroutine cheapest

