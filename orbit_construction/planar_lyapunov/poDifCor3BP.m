 function [x0po,t1] = poDifCor3BP(x0,mu) ;

%        [x0po,t1] = poDifCor3BP(x0,mu) ;
%
% This is the differential correction routine to create a periodic
% orbit (po) about an equilibrium point. It keeps the initial  x
% value constant and varies the y-velocity value.
%
% output: xOpo = initial state on the po (on the negative x-axis) 
%	  t1   = half-period of po
%
% input:  x0   = first guess of initial state on the po 
%         mu   = mass parameter of system
%
%-----------------------------------------------------------------------
%
% Shane Ross (revised 2009-July-07)

global param	% global parameter param is used by integrator

param = mu ;

% tolerances for integration and perpendicular crossing of x-axis

MAXdxdot1 = 1.e-10 ; RelTol = 3.e-14 ; AbsTol = 1.e-16; % high accuracy
% MAXdxdot1 = 1.e-8  ; RelTol = 3.e-10 ; AbsTol = 1.e-12; % low  accuracy
% MAXdxdot1 = 1.e-4  ; RelTol = 3.e-06 ; AbsTol = 1.e-09; % lowest accuracy

MAXattempt = 25;     	% maximum number of attempts before error is declared

dxdot1 	   = 1;   	% to start while loop
attempt    = 0;		% begin counting number of attempts

while abs(dxdot1) > MAXdxdot1
	if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
		break
	end

        % Find first x-axis crossing

        MODEL = 'pcr3bp';
        TSPAN = [0 10] ;  % allow sufficient time for an x-axis crossing
        % high accuracy and Events ON for differential correction
        OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events','on'); 
    	[tt,xx,t1,xx1,i1] = ode113(MODEL,TSPAN,x0,OPTIONS) ;

	x1     = xx1(1) ; 
        y1     = xx1(2) ; 
	dxdot1 = xx1(3) ; 
        ydot1  = xx1(4) ; 

        % Compute the state transition matrix from the initial state to
	% the final state at the next x-axis crossing
      
        % Events option not necessary anymore
        OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 

	[x,t,phi_t1,PHI] = stateTransitionMatrix3BP(x0,t1,mu,OPTIONS) ;

	attempt=attempt+1 ;
	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
	disp(ATTEMPT) ;

	show = 1; % to plot successive orbit  (default=0)
%         show = 0 ; % set to 1 to plot successive approximations
	if show==1,
          plot(x(:,1),x(:,2),'.-'); 
          hold on;
          m = length(x) ;
          plot(x(1,1),x(1,2),'b*');
          plot(x(m,1),x(1,2),'bo');
	  %axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
	  %axis equal
          pause(0.01) ;
        end

   % -------------------------------------
   % this part can be replaced by user for their particular problem

	mu1 = 1-mu; % mass of larger  primary (nearest origin on left)
	mu2 =   mu; % mass of smaller primary (furthest from origin on right)

	r3= ( (x1+mu2)^2 + y1^2 )^1.5;     % r: distance to m1, LARGER MASS
	R3= ( (x1-mu1)^2 + y1^2 )^1.5;     % R: distance to m2, smaller mass

        Ux1 = - x1 + mu1*(x1+mu2)/r3 + mu2*(x1-mu1)/R3 ; % U_x

	xdotdot1 = 2*ydot1 - Ux1 ;  % computing x acceleration

   % -------------------------------------

        dydot0(attempt) = ...
 	    (1/(phi_t1(3,4) - phi_t1(2,4)*(1/ydot1)*xdotdot1)) * dxdot1 ;

        % I added this cheat for large values of the parameter 
        % based on experience with convergence problems
        if mu > 1.e-3,
          DAMP = 1-0.5^attempt ;
        else 
          DAMP = 1 ;
        end

	x0(4) = x0(4) - DAMP*dydot0(attempt) ;
end
x0po=x0 ;
