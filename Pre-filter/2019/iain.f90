program main

    use omp_lib

    implicit none

    real, dimension(4) :: a
    real :: b

    integer :: i

    a(1) = 1.0
    a(2) = 2.0
    a(3) = 3.0
    a(4) = 4.0

    call OMP_SET_NUM_THREADS(2)

    !$OMP PARALLEL DO, PRIVATE(i), REDUCTION(+:b)
    do i = 1,4

        b = b + a(i)

    end do
    !$OMP END PARALEL DO

    WRITE(*,*) "b is", b

end program main
