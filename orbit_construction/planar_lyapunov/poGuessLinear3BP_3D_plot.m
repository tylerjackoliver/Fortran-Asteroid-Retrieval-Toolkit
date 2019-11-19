function poGuessLinear3BP_3D_plot(param,eqNum,Ax,Az,phi,psi) ;

%        poGuessLinear3BP_3D_plot(param,eqNum,Ax,Az,phi,psi) ;
%
% Uses a small displacement from the equilibrium point (in a direction 
% on the collinear point's center manifold) as a first guess for a planar 
% periodic orbit (called a Lyapunov orbit in th rest. three-body problem).
%
% The initial condition and period are to be used as a first guess for
% a differential correction routine.
%
% input:
% param = parameter of system
% eqNum = the number of the equilibrium point of interest
% Ax    = nondim. x-amplitude of periodic orbit 
% Az    = nondim. z-amplitude of periodic orbit 
%
%----------------------------------------------------------------------------
%
% Shane Ross (revised 7.13.04)

mu = param ;	% mass parameter

%eqPos = eqPointLoc3BP(mu,eqNum) ;  % position space location of equil. point
%ep = [eqPos 0 0];		   % phase space location of equil. point

% Get the eigenvalues and eigenvectors of Jacobian of ODEs at equil. point
%[Es,Eu,Ec,Vs,Vu,Vc]=eqPointEig3BP(ep,mu) ;

g = gamma3BP(mu,eqNum) ;   % gamma

c2 = g^(-3)*(mu + (-1)^2*(1-mu)*(g/(1-g))^(2+1)) ;

lam = sqrt( 0.5*(     c2 + sqrt(9*c2^2 - 8*c2) ) ) ;
wp  = sqrt( 0.5*( 2 - c2 + sqrt(9*c2^2 - 8*c2) ) ) 
wv  = c2 ;

kap = (wp^2 + 1 + 2*c2)/(2*wp) ;

i = 0 ;
for t = 0:0.01:pi,
    i = i + 1 ;

    x(i) =    -Ax * cos(wp*t + phi) ;
    y(i) = kap*Ax * sin(wp*t + phi) ;
    z(i) =     Az * sin(wv*t + psi) ;

end


subplot(2,3,4)
plot(x,y,'k') ;
axis equal
AXIS=axis;
XAXIS=1.2*AXIS(1:2);
YAXIS=1.2*AXIS(3:4);
axis([YAXIS YAXIS]);
xlabel('x');
ylabel('y');

subplot(2,3,5)
plot(x,z,'k') ;
axis equal
AXIS=axis;
ZAXIS=1.2*AXIS(3:4);
axis([YAXIS YAXIS]);
xlabel('x');
ylabel('z');

subplot(2,3,6)
plot(y,z,'k') ;
axis equal
axis([YAXIS YAXIS]);
xlabel('y');
ylabel('z');

subplot(2,3,2)
plot3(x,y,z,'k');
axis equal
axis([YAXIS YAXIS YAXIS]);
xlabel('x');
ylabel('y');
zlabel('z');

