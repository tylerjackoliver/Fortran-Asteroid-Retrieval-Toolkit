function state_out = rotator(state_in, epoch, tt)

    GMSun = 1.32712440018e11; %km^3/s^2

    elts = cspice_oscelt(state_in, 0., GMSun); % States are defined at J2000!
    
    % Apply relations from Sanchez et. al. The following code
    % block assumes that the time t is inputted as ephemeris
    % seconds past J2000 (and thus that it is linked innately to
    % the epoch of the states used in this work.
   
    if (abs(elts(4)) < 1.0e-06)

        elts(4) = 0;
        elts(5) = elts(5) + (epoch+tt)*(2.*pi)/(365.25*86400.);                                                                   ! Add on phasing equal to the transfer time                                           

    else

        elts(4) = elts(4) + (epoch+tt)*(2.d0 * pi)/(365.25d0*86400.d0);                                                                  ! Add on phasing equal to the transfer time
        elts(5) = 0;

    end
    
    state_out = cspice_conics(elts, epoch+tt);

end