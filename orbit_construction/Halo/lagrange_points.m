function [L] = lagrange_points(mu)
%
% File: lagrange_points.m
% Author: Nick Truesdale
% Date: 12/10/2012
%
% Description: This function numerically solves equations from the CRTBP to
% find the five Lagrange points in the synodic system. They are returned in
% L as a 5x2 array of (x,y) points.

% Collinear equation
f = @(x, a, b) x.^2.*( (1-b) + 3*x + 3*x.^2 + x.^3 ) - ...
               mu*(a + 2*a*x + (1 + a - b)*x.^2 + 2*x.^3 + x.^4 );

% Options
opt = optimset('display', 'off', 'tolfun', 1e-14, 'tolx', 1e-14, ...
               'maxiter', 1000);

% Preallocate L
L = zeros(5,2);
           
% Solve
L(1) = fsolve(@(x) f(x, -1,  1),  1, opt) + 1 - mu;
L(2) = fsolve(@(x) f(x,  1,  1),  1, opt) + 1 - mu;
L(3) = -(fsolve(@(x) f(x, -1, -1), 1, opt) + 1 - mu);

% Points 4 and 5
L(4,:) = [0.5 - mu, sqrt(3)/2];
L(5,:) = [0.5 - mu, -sqrt(3)/2];

end % function