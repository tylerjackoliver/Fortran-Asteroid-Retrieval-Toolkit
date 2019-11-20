
subroutine cheaper(a_b, e_b, inb, rp_t, e_t, in_t, M_t, cheap)

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

    ! First case: perihelion modifying aphelion

    ! a_int = (rp_b+rt)/2

    ! dV1 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    ! dV2 = sqrt(mu_s*(2/rt-1/a_t))-sqrt(mu_s*(2/rt-1/a_int))

    ! dVi1 = 2*sqrt((mu_s/a_b)*(ra_b/rp_b))*sin(abs(in_b-in_t)/2)

    ! if (r+)


    ! a_t = rp_t/(1-e_t) 
    ! ra_t = a_t*(1+e_t) 

    ! There are 16 cases of different transfer possibilities. 

    !! Cases 1-4: target rt, perform periapsis changing apoapsis first

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


    ! dVi8 = 2*sqrt((mu_s/a_t)*ra_t/rp_t)*sin(abs(in_b-in_t)/2) 

    v5 = sqrt(dVi5**2+dV3**2)+sqrt(dV4**2) 
    v6 = sqrt(dVi6**2+dV3**2)+sqrt(dV4**2)
    v7 = sqrt(dV3**2) + sqrt(dV4**2+dVi7**2) 

    ! !! Cases 9-12: target ra_t, perform periapsis changing apoapsis first

    ! a_int = .5*(ra_t+rp_b) 

    ! dV5 = sqrt(mu_s*(2/rp_b-1/a_int))-sqrt(mu_s*(2/rp_b-1/a_b)) 
    ! dV6 = sqrt(mu_s*(2/ra_t-1/a_t)) - sqrt(mu_s*(2/ra_t-1/a_int)) 

    ! dVi9 = 2*sqrt((mu_s/a_b)*(ra_b/rp_b))*sin(abs(in_b-in_t)/2) 

    ! if (ra_t > rp_b) then

    !     dVi10 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
    !     dVi11 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
        
    ! else

    !     dVi10 = 2*sqrt((mu_s/a_int)*rp_b/ra_t)*sin(abs(in_b-in_t)/2) 
    !     dVi11 = 2*sqrt((mu_s/a_int)*ra_t/rp_b)*sin(abs(in_b-in_t)/2) 
        
    ! end if

    ! dVi12 = 2*sqrt((mu_s/a_t)*(rp_t)/ra_t)*sin(abs(in_b-in_t)/2) 

    ! v9 = sqrt(dVi9**2+dV5**2)+sqrt(dV6**2) 
    ! v10 = sqrt(dV5**2+dVi10**2)+sqrt(dV6**2) 
    ! v11 = sqrt(dV5**2) + sqrt(dV6**2+dVi11**2) 
    ! v12 = sqrt(dV5**2) + sqrt(dV6**2+dVi12**2) 

    ! !! Cases 13-16: target ra_t, perform ap changing peri first

    ! a_int = .5*(ra_b+ra_t) 

    ! dV7 = sqrt(mu_s*(2/ra_b-1/a_int))-sqrt(mu_s*(2/ra_b-1/a_b)) 
    ! dV8 = sqrt(mu_s*(2/ra_t-1/a_t))-sqrt(mu_s*(2/ra_t-1/a_int)) 

    ! dVi13 = 2*sqrt(mu_s/a_b*(rp_b/ra_b))*sin(abs(in_b-in_t)) 

    ! if (ra_b > ra_t) then

    !     dVi14 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 
    !     dVi15 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 

    ! else

    !     dVi14 = 2*sqrt(mu_s/a_int*ra_b/ra_t)*sin(abs(in_b-in_t)) 
    !     dVi15 = 2*sqrt(mu_s/a_int*ra_t/ra_b)*sin(abs(in_b-in_t)) 

    ! end if

    ! dVi16 = 2*sqrt(mu_s/a_int*(ra_t/rp_t))*sin(abs(in_b-in_t)) 

    ! v13 = sqrt(dVi13**2+dV7**2)+sqrt(dV8**2) 
    ! v14 = sqrt(dV7**2+dVi14**2)+sqrt(dV8**2) 
    ! v15 = sqrt(dV7**2) + sqrt(dV8**2+dVi15**2) 
    ! v16 = sqrt(dV7**2) + sqrt(dV8**2+dVi16**2) 

    cheap = min(v1, v2, v3, v5, v6, v7)

end subroutine cheaper
