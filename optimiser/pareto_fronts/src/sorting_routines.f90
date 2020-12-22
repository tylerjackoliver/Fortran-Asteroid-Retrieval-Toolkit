! @brief Contains subroutines from Numerical Recipes in Fortran for sorting arrays rapidly

module sorting_routines

    implicit none

    contains

    ! @brief Sorts a 1D array in Fortran from negative to positive
    ! @param[in] n Size of array
    ! @param[inout] Array Array to be sorted; on exit, it is the sorted array
    ! @param[inout] Index Permutation matrix for the array sorting
    subroutine sort_array(n, Array, Index)
            
        implicit none
        
        integer,            intent(in)  :: n
        double precision,   intent(in)  :: Array(n)
        integer,            intent(out) :: Index(n)
        integer, parameter              :: nn=15, nstack=1000
        integer                         :: k, i, j, indext, jstack, l, r
        integer                         :: istack(nstack)
        double precision                :: a
        
        do j = 1,n
            Index(j) = j
        end do
        jstack=0
        l=1
        r=n
        do
            if (r-l < nn) then
                do j=l+1,r
                    indext=Index(j)
                    a=Array(indext)
                    do i=j-1,l,-1
                        if (Array(Index(i)) <= a) exit
                        Index(i+1)=Index(i)
                    end do
                    Index(i+1)=indext
                end do
                if (jstack == 0) return
                r=istack(jstack)
                l=istack(jstack-1)
                jstack=jstack-2
            else
                k=(l+r)/2
                call swap(Index(k),Index(l+1))
                call exchangeIndex(Index(l),Index(r))
                call exchangeIndex(Index(l+1),Index(r))
                call exchangeIndex(Index(l),Index(l+1))
                i=l+1
                j=r
                indext=Index(l+1)
                a=Array(indext)
                do
                    do
                        i=i+1
                        if (Array(Index(i)) >= a) exit
                    end do
                    do
                        j=j-1
                        if (Array(Index(j)) <= a) exit
                    end do
                    if (j < i) exit
                    call swap(Index(i),Index(j))
                end do
                Index(l+1)=Index(j)
                Index(j)=indext
                jstack=jstack+2
                if (jstack > nstack) then
                    write(*,*) 'NSTACK too small in indexArrayReal()'   ! xxx
                    error stop
                end if
                if (r-i+1 >= j-l) then
                    istack(jstack)=r
                    istack(jstack-1)=i
                    r=j-1
                else
                    istack(jstack)=j-1
                    istack(jstack-1)=l
                    l=i
                end if
            end if
        end do
    contains
        subroutine exchangeIndex(i,j)
            integer, intent(inout) :: i,j
            integer                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
        pure elemental subroutine swap(a,b)
            implicit none
            integer, intent(inout) :: a,b
            integer :: dum
            dum=a
            a=b
            b=dum
        end subroutine swap
    end subroutine sort_array
end module sorting_routines