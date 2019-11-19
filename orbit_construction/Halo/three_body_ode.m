function [out1, isterminal, direction] = three_body_ode(t, y, flag)
%
% File: three_body_ode.m
% Author: Nick Truesdale
% Date: 12/10/2012
%
% Description: This function gives the equations of motion for the CRTBP
% in a rotating synodic frame. The state Y contains three-space position 
% and velocity coordinates of the third body. MU is the mass parameter.

mu = 3.040402499929356e-06;

% Vectors between primaries and satellite
r1 = y(1:3);
r1(1) = r1(1) + mu;

r2 = r1;
r2(1) = r2(1) - 1;

% Vector norms
R1 = norm(r1)^3;
R2 = norm(r2)^3;

% Acceleration derivatives
yprime(4:6,1) = (mu - 1)*r1/R1 - mu*r2/R2;
yprime(4:6,1) = yprime(4:6,1) + [2*y(5) + y(1);  
                                -2*y(4) + y(2);
                                            0];
% Velocity derivatives
yprime(1:3,1) = y(4:6);

out1 = yprime;

switch flag
    
    case 'events'
        
        isterminal = 1;
        out1 = atan2(y(2),y(1)) - pi/8;
        direction = 0;
        
end

end