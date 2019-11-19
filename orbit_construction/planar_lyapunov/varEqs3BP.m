function PHIdot = varEqs3BP(t,PHI) ;

%        PHIdot = varEqs3BP(t,PHI) ;
%
% This here is a preliminary state transition, PHI(t,t0),
% matrix equation attempt for the planar CR3BP, based on...
%
%        d PHI(t, t0)
%        ------------ =  Df(t) * PHI(t, t0)
%             dt
%
%-----------------------------------------------------------
% CR3BP CONVENTION:
%                 L4
%
%
%    L3-----M1-------L1---M2---L2         M1=1-mu, M2=mu
%
%
%                 L5
%
% Shane Ross (revised 2.19.04)

global param

% mu = mass paramater

mu = param ; 

mu1 = 1-mu ;
mu2 =   mu ;

x(1:4) = PHI(17:20);
phi  = reshape(PHI(1:16),4,4);

r2= (x(1)+mu2)^2 + x(2)^2;      % r: distance to m1, LARGER MASS
R2= (x(1)-mu1)^2 + x(2)^2;      % R: distance to m2, smaller mass

r3= r2^1.5;
r5= r2^2.5;
R3= R2^1.5;
R5= R2^2.5;

% The following are the two partial derivatives of the
% effective potential U(x,y)

Ux = - x(1) + mu1*(x(1)+mu2)/r3 + mu2*(x(1)-mu1)/R3 ;
Uy = - x(2) + mu1* x(2)     /r3 + mu2* x(2)     /R3 ;

% The following are three double partial derivatives of the
% effective potential U(x,y)

Uxx = -1+(mu1/r3)*(1-(3*(x(1)+mu2)^2/r2))+(mu2/R3)*(1-(3*(x(1)-mu1)^2/R2)) ;
Uyy = -1+(mu1/r3)*(1-(3* x(2)     ^2/r2))+(mu2/R3)*(1-(3* x(2)     ^2/R2)) ;
Uxy =   -(mu1/r5)*    3* x(2)*(x(1)+mu2) -(mu2/R5)*    3* x(2)*(x(1)-mu1)  ;

% The following is the Jacobian matrix

Df    =[  0     0    1    0 ;
          0     0    0    1 ;
        -Uxx  -Uxy   0    2 ;
        -Uxy  -Uyy  -2    0 ];

phidot = Df * phi; % variational equation

PHIdot        = zeros(20,1);
PHIdot(1:16)  = reshape(phidot,16,1);
PHIdot(17)    = x(3);
PHIdot(18)    = x(4);
PHIdot(19)    = 2*x(4) - Ux ;
PHIdot(20)    =-2*x(3) - Uy ;
