function [xf, tf, success] = differentialCorrector(x0, mu)

    tolerance = 1e-012;
    OPTIONS_events = odeset('RelTol', tolerance, 'AbsTol', tolerance, ...
                            'Events', @findZ);
                 
    % Set iteration counters and maximums
    
    max_iterations = 200;
    iteration_count = 0;
    
    success = false;                                                        % Success flag
    
    while (iteration_count < max_iterations)
        
       y0 = [x0; reshape(eye(6), 36, [])];                                  % Add STM to state vector
       
       [~, ~, te, ye, ~] = ode45(@cr3bp, [0 10], y0, OPTIONS_events, mu);   % Integrate to z-axis crossing
       
       % Periodicity check - proceed no further if all is OK
       
       if (abs(ye(2)) < 1e-07 && abs(ye(4)) < 1e-07)
          
           success = true;
           break;
           
       end
       
       % Extract end conditions
       
       x = ye(1);
       y = ye(2);
       z = ye(3);
       xdot = ye(4);
       ydot = ye(5);
       zdot = ye(6);
       
       % Compute radial distances
       
       r1 = sqrt((x+mu)^2 + y^2 + z^2);	% Distance to primary
       r2 = sqrt((x-1+mu)^2 + y^2+z^2);	% Distance to secondary
       
       % Compute x-acceleration
       
       coeff1 = (1-mu)/r1^3;                                                % Define coefficients to make our life easier
       coeff2 = mu/r2^3;                                                    % ...
       
       xddot = 2 * ydot + x - coeff1 * (x+mu) - coeff2 * (x-1+mu);
           
       % Extract the state transition matrix
       
       STM = reshape(ye(7:end), 6, []);
       
       % Compute corrections
       
       corr = ([STM(2,3) STM(2,5); STM(4,3) STM(4,5)]- ...
              (1/zdot)*[ydot; xddot]*[STM(3,3) STM(3,5)])\[-y; -xdot];
       
       % Update velocities
       
       x0(3) = x0(3) + corr(1);
       x0(5) = x0(5) + corr(2);
       
       % Iterate counter
       
       iteration_count = iteration_count + 1;
               
    end
    
    xf = x0;
    tf = 4*te;

end

%% Dependency subroutines

%
%	Generates EOM, using variational equations & STM
%

function xdot = cr3bp(~, x, mu)

	dimensions = 6;                                                         % Set up number of dimensions of phase space
	x0 = x(1:dimensions)';                                                  % Extract trajectory from propagated x vector

	phi = reshape(x(dimensions+1:end), ...
     dimensions, []);                                                       % Extract STM

	jacobian = computeJacobian(x0, mu);                                     % Compute the jacobian derivative matrix

	stmdot = reshape((jacobian*phi), [], 1);                                % Grab the derivative of the STM
	
    dx = zeros(6, 1);                                                       % Initialise state derivative matrix
	
    r1 = sqrt((x0(1)+mu)^2 + x0(2)^2 + x0(3)^2);                            % Distance to primary
	r2 = sqrt((x0(1)-1+mu)^2 + x0(2)^2+x0(3)^2);                            % Distance to secondary
	
    dx(1:3) = x0(4:6);                                                      % d/dt form first three elements of derivative matrix
	
    coeff1 = (1-mu)/r1^3;                                                   % Define coefficients to make our life easier
	coeff2 = mu/r2^3;				
	
    dx(4) = 2*x0(5)+x0(1)-coeff1*(x0(1)+mu)-...
	coeff2*(x0(1)-1+mu);                                                    % x-accel using potential function
	
    dx(5) = -2*x0(4)+x0(2)-coeff1*x0(2)...
	-coeff2*x0(2);                                                          % y-accel
	
    dx(6) = -coeff1*x0(3)-coeff2*x0(3);                                     % z-accel
	
    xdot = [dx; stmdot];                                                    % Concatenate derivatives

end

%
%	Generates jacobian matrix for computing STM derivative
%	
%	Contains variational equations for three-dimensions; see Parker (2014)
%

function jacobian = computeJacobian(x0, mu)
	
    r1 = sqrt((x0(1)+mu)^2 + x0(2)^2 + x0(3)^2);                            % Distance to primary
	r2 = sqrt((x0(1)-1+mu)^2 + x0(2)^2+x0(3)^2);                            % Distance to secondary
	
    j = [0 1 0; -1 0 0; 0 0 0];                                             % First component of jacobian
	
    % Partial derivatives - start with unique terms...
	
    Uxx = 1-(1-mu)/r1^3-mu/r2^3+...
          3*(1-mu)*(x0(1)+mu)^2/r1^5+...
          3*mu*(x0(1)-1+mu)^2/r2^5;
	
    Uyy = 1-(1-mu)/r1^3-mu/r2^3+...
          3*(1-mu)*x0(2)^2/r1^5+3*mu*x0(2)^2/r2^5;
	
    Uzz = (1-mu)/r1^3-mu/r2^3+3*(1-mu)*x0(3)^2/r1^5+...
          3*mu*x0(3)^2/r2^5;
	
    Uxy = 3*(1-mu)*(x0(1)+mu)*x0(2)/r1^5+...
          3*mu*(x0(1)-1+mu)*x0(2)/r2^5;
	
    Uxz = 3*(1-mu)*(x0(1)+mu)*x0(3)/r1^5+...
          3*mu*(x0(1)-1+mu)*x0(3)/r2^5;
	
    Uyz = 3*(1-mu)*x0(2)*x0(3)/r1^5+...
          3*mu*x0(2)*x0(3)/r2^5;
	
    % And now the repeated terms...
	
    Uyx = Uxy;
	Uzx = Uxz;
	Uzy = Uyz;
	
    uMatrix = [Uxx Uxy Uxz; ...
               Uyx Uyy Uyz; ...
               Uzx Uzy Uzz];                                                % Set up u-matrix
	jacobian = [zeros(3), eye(3); uMatrix, 2*j];                            % And construct the jacobian

end

function [value, isterminal, direction] = findZ(~, x, ~)

    value = x(3);
    isterminal = 1;
    direction = -1;

end
