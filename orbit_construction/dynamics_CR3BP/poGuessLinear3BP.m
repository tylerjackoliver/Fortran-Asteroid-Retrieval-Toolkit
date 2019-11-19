function [x0poGuess,TGuess] = poGuessLinear3BP(param,eqNum,Ax) ;

%        [x0poGuess,TGuess] = poGuessLinear3BP(param,eqNum,Ax) ;
%
% Uses a small displacement from the equilibrium point (in a direction 
% on the collinear point's center manifold) as a first guess for a planar 
% periodic orbit (called a Lyapunov orbit in th rest. three-body problem).
%
% The initial condition and period are to be used as a first guess for
% a differential correction routine.
%
% output:
% x0poGuess = initial state on the periodic orbit (first guess)
%           = [ x 0  0 yvel]  , (i.e., perp. to x-axis and in the plane)
% TGuess    = period of periodic orbit (first guess)
%
% input:
% param = parameter of system
% eqNum = the number of the equilibrium point of interest
% Ax    = nondim. x-amplitude of periodic orbit 
%
%----------------------------------------------------------------------------
%
% Shane Ross (revised 2.17.04)

mu = param ;	% mass parameter

eqPos = eqPointLoc3BP(mu,eqNum) ;  % position space location of equil. point
ep = [eqPos 0 0];		   % phase space location of equil. point

% Get the eigenvalues and eigenvectors of Jacobian of ODEs at equil. point
[Es,Eu,Ec,Vs,Vu,Vc]=eqPointEig3BP(ep,mu) ;

l = abs(imag(Ec(1))) ;  % lambda

g = gamma3BP(mu,eqNum) ;  % gamma

c2 = g^(-3)*(mu + (-1)^2*(1-mu)*(g/(1-g))^(2+1)) ;

k = (l^2 + 1 + 2*c2)/(2*l) ;

x0poGuess	= zeros(4,1) ;
x0poGuess(1)	= ep(1) - Ax ;
x0poGuess(4)	= ep(4) + Ax*k*l ;

TGuess = 2*pi/l ;
