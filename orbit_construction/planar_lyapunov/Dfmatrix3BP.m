function Df = Dfmatrix3BP(x,mu)

%        Df = Dfmatrix3BP(x,mu);
%
% Gives matrix Df(x) (i.e., the matrix of derivatives of f, where xdot=f(x) )
% for a 4-dimensional point x in the planar CR3BP's phase space 
%
%-----------------------------------------------------------------------
% CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
% with the LARGER MASS, M1 to the left of the origin at (-mu,0)
% and the smaller mass, M2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%       (rotating coords)
%
%                 L4
%
% -L3------M1--+-----L1--M2--L2-
%
%                 L5
%
% Shane Ross (revised 2.19.04)

% mu = mass paramater

mu1 = 1-mu ;
mu2 =   mu ;

r2= (x(1)+mu2)^2 + x(2)^2;      % r: distance to m1, LARGER MASS
R2= (x(1)-mu1)^2 + x(2)^2;      % R: distance to m2, smaller mass

r3= r2^1.5;
r5= r2^2.5;
R3= R2^1.5;
R5= R2^2.5;

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
