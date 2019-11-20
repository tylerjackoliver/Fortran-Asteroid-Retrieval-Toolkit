% Rotates conditions from the CR3BP synodic frame into the global frame

function stateOut = invRotationMatrix_final(state, t)

mu = 3.03e-06;

r = state(1:3);
r(1) = r(1) + mu;
rdot = state(4:6);

theta0 = 100.3762*pi/180;
total_angle = theta0 + (t-51544.5)/365.25*2*pi; % t is linked to the angle, so just use t

t_ir = [cos(total_angle) -sin(total_angle) 0; ...
        sin(total_angle) cos(total_angle) 0; ...
        0 0 1];

t_irdot = [-sin(total_angle) -cos(total_angle) 0;...
           cos(total_angle) -sin(total_angle) 0; ...
           0 0 0];
       
rPrime = t_ir*r';
vPrime = t_ir * rdot' + t_irdot * r';

stateOut = [rPrime; vPrime]';
       
end