function [T] = richardson (mu,Lpt,Azlp,n);
% [T] = richardson (mu,Lpt,Azlp,n);
%
% Gives the state (position,velocity) for a 3D periodic halo orbit 
% centered on the specified collinear Lagrange point as a function
% of an angular variable 
% [Uses Richardson's 3rd order model for analytically constructing a 3D 
%  periodic halo orbit about the points L1, L2, or L3]
% 
% output: xx = (r,v) state desired halo orbit (in 3D CR3BP nondim. units)
%	  th = angular variable from 0 to 2*pi (about 100 points)
%
% input: mu   = mass parameter of system [ mu = m2/(m1+m2) ]
% 	 Lpt  = 1,2, or 3, the number of the specified collinear Lagrange point
% 	 Azlp = out-of-plane (or z-amplitude) of the desired halo 
%	        [ in Lpt-primary distances (see gammaL.m) ] 
% 	 n    = +1 is northern halo (z>0, Class I), 
% 	      = -1 is southern halo (z<0, Class II)
%
%----------------------------------------------------------------------------
% CR3BP with the LARGER MASS, m1 to the left of the origin at (-mu,0)
% and m2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%                L4
%
% -L3----m1--+-------L1--m2--L2-
%
%                L5
%
% Shane Ross (revised 7.13.04)

Az = Azlp;

gamma = gamma3BP(mu,Lpt) ;   % gamma

if     Lpt == 1, won =+1; primary = 1-mu;
elseif Lpt == 2, won =-1; primary = 1-mu;
elseif Lpt == 3, won =+1; primary =  -mu; end

if Lpt == 3
    for N = 2:4
	c(N) = (1/gamma^3)*( 1-mu + (-primary*gamma^(N+1))/((1+gamma)^(N+1)) );
    end
else
    for N = 2:4
	c(N) = (1/gamma^3)*( (won^N)*mu + ...
	       ((-1)^N)*((primary)*gamma^(N+1))/((1+(-won)*gamma)^(N+1)) );
    end
end

polylambda = [ 1 0 (c(2)-2) 0 -(c(2)-1)*(1+2*c(2)) ];
lambda = sort(roots(polylambda)); % lambda = frequency of orbit
if Lpt == 3, lambda = abs(lambda(3)) ;
else,	     lambda = abs(lambda(1)) ; 
end

k = 2*lambda/(lambda^2 + 1 - c(2));
del = lambda^2 - c(2);

d1 = ((3*lambda^2)/k)*(k*(6*lambda^2 - 1) - 2*lambda);
d2 = ((8*lambda^2)/k)*(k*(11*lambda^2 -1) - 2*lambda);

a21 = (3*c(3)*(k^2 - 2))/(4*(1 + 2*c(2)));
a22 = 3*c(3)/(4*(1 + 2*c(2)));
a23 = -(3*c(3)*lambda/(4*k*d1))*( 3*(k^3)*lambda - 6*k*(k-lambda) + 4);
a24 = -(3*c(3)*lambda/(4*k*d1))*( 2 + 3*k*lambda );

b21 = -(3*c(3)*lambda/(2*d1))*(3*k*lambda - 4);
b22 = 3*c(3)*lambda/d1;
d21 = -c(3)/(2*lambda^2);

a31 = -(9*lambda/(4*d2))*(4*c(3)*(k*a23 - b21) + k*c(4)*(4 + k^2)) + ((9*lambda^2 + 1 -c(2))/(2*d2))*(3*c(3)*(2*a23 - k*b21) + c(4)*(2 + 3*k^2));
a32 = -(1/d2)*( (9*lambda/4)*(4*c(3)*(k*a24 - b22) + k*c(4)) + 1.5*(9*lambda^2 + 1 - c(2))*( c(3)*(k*b22 + d21 - 2*a24) - c(4)) );

b31 = (.375/d2)*( 8*lambda*(3*c(3)*(k*b21 - 2*a23) - c(4)*(2 + 3*k^2)) + (9*lambda^2 + 1 + 2*c(2))*(4*c(3)*(k*a23 - b21) + k*c(4)*(4 + k^2)) );
b32 = (1/d2)*( 9*lambda*(c(3)*(k*b22 + d21 - 2*a24) - c(4)) + .375*(9*lambda^2 + 1 + 2*c(2))*(4*c(3)*(k*a24 - b22) + k*c(4)) );

d31 = (3/(64*lambda^2))*(4*c(3)*a24 + c(4));
d32 = (3/(64*lambda^2))*(4*c(3)*(a23 - d21) + c(4)*(4 + k^2));

s1 = (1/(2*lambda*(lambda*(1+k^2) - 2*k)))*( 1.5*c(3)*(2*a21*(k^2 - 2)-a23*(k^2 + 2) - 2*k*b21) - .375*c(4)*(3*k^4 - 8*k^2 + 8) );
s2 = (1/(2*lambda*(lambda*(1+k^2) - 2*k)))*( 1.5*c(3)*(2*a22*(k^2 - 2)+a24*(k^2 + 2) + 2*k*b22 + 5*d21) + .375*c(4)*(12 - k^2) );

a1 = -1.5*c(3)*(2*a21+ a23 + 5*d21) - .375*c(4)*(12-k^2);
a2 = 1.5*c(3)*(a24-2*a22) + 1.125*c(4);

l1 = a1 + 2*(lambda^2)*s1;
l2 = a2 + 2*(lambda^2)*s2;

tau1 = 0:0.01:2*pi+0.01 ;
deltan =-n;

Ax = sqrt( (-del - l2*Az^2)/l1 );

x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau1) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);
y = k*Ax*sin(tau1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau1) + (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau1);
z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau1);
xdot = lambda*Ax*sin(tau1) - 2*lambda*(a23*Ax^2 - a24*Az^2)*sin(2*tau1) - 3*lambda*(a31*Ax^3 - a32*Ax*Az^2)*sin(3*tau1);
ydot = lambda*(k*Ax*cos(tau1) + 2*(b21*Ax^2 - b22*Az^2)*cos(2*tau1) + 3*(b31*Ax^3 - b32*Ax*Az^2)*cos(3*tau1));
zdot = - lambda*deltan*Az*sin(tau1) - 2*lambda*deltan*d21*Ax*Az*sin(2*tau1) - 3*lambda*deltan*(d32*Az*Ax^2 - d31*Az^3)*sin(3*tau1);

r0 = gamma*[ (primary + gamma*(-won +x(:)))/gamma -y(:) z(:) ];
v0 = gamma*[ xdot(:) ydot(:) zdot(:) ];
xx= [r0 v0];

th = tau1(:) ;

gam = 1.497531218885853e+06 ; % scale factor for Sun-Earth L1
