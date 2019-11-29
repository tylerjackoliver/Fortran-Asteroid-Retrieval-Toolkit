function stateOut = globalTocr3bp(state, t)

mu = 3.0003e-06;

% Convert t in seconds to t in days past J2000 (from ephemeris time, so t
% is linked to j2000

t = t / (365.25 * 86400.);

r = state(1:3);
rdot = state(4:6);
theta0 = 100.3762*pi/180;
total_angle = theta0 + (t-51544.5)/365.25*2*pi; % t is linked to the angle, so just use t

t_ir = [cos(total_angle) sin(total_angle) 0; ...
        -sin(total_angle) cos(total_angle) 0; ...
        0 0 1];

t_irdot = [-sin(total_angle) cos(total_angle) 0;...
           -cos(total_angle) -sin(total_angle) 0; ...
           0 0 0];
       
rPrime = t_ir*r;
vPrime = t_ir * rdot + t_irdot * r;

stateOut = [rPrime; vPrime];

stateOut(1) = stateOut(1) - mu;
       
end
