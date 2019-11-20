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
    integer :: desiredbody

    real :: a_b, e_b, inb, rp_t, e_t, in_t, start, finish, M_t
    real :: mindV, au,  prevCheapest
    real :: t1, t2, t3, t4, t5, date, o, Om, M, p


    real, dimension(1353, 3) :: bodyData
    real, dimension(3046544,4) :: targetData
    real, dimension(1,5) :: minConds

    character(len=260) :: bodystring, filepath
    character(len=8) :: fmt = '(I7.7)'

    candidatecount = 0
    itercount = 0
    au =  1.495978707e8
    filepath = 'corrected_oe/'
    prevCheapest = 10000

    ! Print welcome message

    write(*,*) "============================================================"
    write(*,*) " "
    write(*,*) "OPTIMAL MANIFOLD CALCULATOR"
    write(*,*) " "
    write(*,*) " "
    write(*,*) "============================================================"

    ! Initialise OMP and delete current results.txt

    call OMP_SET_NUM_THREADS(nthreads)    
    call file_check()
    ! call getarg(1, desiredbody)

    ! Load data sets

    write(*,*) "Loading data sets..."

    open(unit=3, file="target_csv.csv")    ! Target OEs

    do i = 1,3046544

        read(3,*) targetData(i,1), targetData(i,2), targetData(i,3),t1,t2,targetData(i,4),t4,t5

    end do

    close(3)

    ! Write progress message and begin timing

    write(*,*) "Loading complete. Beginning compute sequence for body #", desiredbody
    start = omp_get_wtime()

        bodyID = 3403148
        write(bodystring, fmt) bodyID

        unitno = OMP_GET_THREAD_NUM() + 30

        open(unit=unitno, file=trim(filepath)//trim(bodystring)//'_OE_corr.csv')
    
        do j = 1,1353
            read(unitno,*) date, bodyData(j,1), bodyData(j,2), bodyData(j,3), o, O, M, p
        end do

        close(unitno)


    	!$OMP PARALLEL DO PRIVATE(I,J,K, bodyid, bodystring, unitno, mindV, bodyData, a_b, e_b, inb, rp_t, e_t, in_t)	
        do j = 272,1353

            a_b = bodyData(j,1)*au
            e_b = bodyData(j,2)
            inb = bodyData(j,3)

            do k = 1,3046544
                
                rp_t = targetData(k,1)*au
                e_t = targetData(k,2)
                in_t = targetData(k,3)
                M_t = targetData(k,4)

                call cheapest(a_b, e_b, inb, rp_t, e_t, in_t, M_t, mindV)

                if (mindV < prevCheapest) then 

                    minConds(1,1) = rp_t
		    minConds(1,2) = e_t
		    minConds(1,3) = in_t
		    minConds(1,4) = M_t
		    minConds(1,5) = mindV

                end if

            end do

        write(*,*) "Completed processing date #",  j

        end do

    end do
    !$OMP END PARALLEL DO

    finish = omp_get_wtime()

    write(*,*) "==================== Scan completed ========================"
    write(*,*) "Minimum dV: ", minConds(1,5)
    write(*,*) "Radius of periapsis: ", minConds(1,1)
    write(*,*) "Eccentricity: ", minConds(1,2)
    write(*,*) "Inclination [rad]: ", in_t
    write(*,*) "Mean anomaly [rad]: ", M_t
    write(*,*) "============================================================"
v

end program main


     subroutine cheapest(a_b, e_b, inb, rp_t, e_t, in_t, M_t, cheap)

         real, intent(in) :: a_b, e_b, inb, rp_t, e_t, in_t, M_t
         real, intent(out) :: cheap
         real :: mu_s, a_t, ra_t, rp_b, ra_b, a_int, in_b
         real :: dV1, dV2, dV3, dV4
         real :: dVi1, dVi2, dVi3, dVi5, dVi6, dVi7, dVi4, dVi8
         real :: v1, v2, v3, v5, v6, v7
         double precision :: pi = 4.*datan(1.d0)
         real :: rt

         mu_s = 1.327124400189e11

         in_b = inb*pi/180

         ra_b = a_b*(1+e_b)    ! if (r+)
         rp_b = a_b*(1-e_b)

         a_t = rp_t/(1-e_t)
         rt = a_t*(1-e_t**2)/(1+e_t*cos(M_t))
         ra_t = a_t*(1+e_t)

         a_int = .5*(rp_b+rt) 

         dV1 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
         dV2 = sqrt(mu_s*(2/rt-1/a_t)) - sqrt(mu_s*(2/rt-1/a_int)) 

         dVi1 = 2*sqrt((mu_s/a_b)*(ra_b/rp_b))*sin(abs(in_b-in_t)/2)

         if (rt>rp_b) then


         dVi2 = 2*sqrt((mu_s/a_int) * (rt/rp_b))*sin(abs(in_b-in_t)/2)
         dVi3 = 2*sqrt((mu_s/a_int) * (rp_b/rt))*sin(abs(in_b-in_t)/2) 

         else

         dVi2 = 2*sqrt((mu_s/a_int)*(rp_b/rt))*sin(abs(in_b-in_t)/2)
         dVi3 = 2*sqrt((mu_s/a_int)*(rt/rp_b))*sin(abs(in_b-in_t)/2) 

         end if

         ! dVi4 = 2*sqrt((mu_s/a_t)*(rp_b/rt))*sin(abs(in_b-in_t)/2) 

         v1 = sqrt(dV1**2+dVi1**2) + sqrt(dV2**2)
         v2 = sqrt(dV1**2+dVi2**2) + sqrt(dV2**2) 
         v3 = sqrt(dV1**2) + sqrt(dVi3**2+dV2**2) 

         ! !! Cases 5:8: target rt, perform apoapsis changing periapsis first

         a_int = .5*(rt+ra_b) 

         dVi5 = 2*sqrt((mu_s/a_b)*(rp_b/ra_b))*sin(abs(in_b-in_t)/2)

         dV3 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b))
         dV4 = sqrt(mu_s*(2/rt-1/a_t))-sqrt(mu_s*(2/rt-1/a_int))

         if (rt > ra_b) then

         dVi6 = 2*sqrt((mu_s/a_int)*(ra_b/rt))*sin(abs(in_b-in_t)/2)
         dVi7 = 2*sqrt((mu_s/a_int)*(rt/ra_b))*sin(abs(in_b-in_t)/2)

         else

         dVi6 = 2*sqrt((mu_s/a_int)*(rt/ra_b))*sin(abs(in_b-in_t)/2) 
         dVi7 = 2*sqrt((mu_s/a_int)*(ra_b/rt))*sin(abs(in_b-in_t)/2)

         end if

         v5 = sqrt(dVi5**2+dV3**2)+sqrt(dV4**2) 
         v6 = sqrt(dVi6**2+dV3**2)+sqrt(dV4**2)
         v7 = sqrt(dV3**2) + sqrt(dV4**2+dVi7**2) 

         cheap = min(v1, v2, v3, v5, v6, v7)

     end subroutine cheapest

subroutine data_index(k, outval1, outval2)

	integer, intent(in) :: k
	integer, intent(out) :: outval1, outval2 


	    if (k <= 2791440) then
	
		outval1 = floor(REAL(k)/360)*360
		outval2 = outval1 + 360

	    else if (k > 2791440 .AND. k <= 2800674) then
	
		outval1 = floor((REAL(k)-2791440)/32)*32+2791440
		outval2 = outval1+32

	    else if (k > 2800674) then

		outval1 = floor((REAL(k)-2800674)/360)*360+2800674
		outval2 = outval1+360
	
	    end if

end subroutine data_index
