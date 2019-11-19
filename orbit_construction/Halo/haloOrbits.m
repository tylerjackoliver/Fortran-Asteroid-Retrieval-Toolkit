% Function to generate corrected three-dimensional halo orbits in the CR3BP

%% Initialisation

mu = 3.003458e-06;                                                          % Mass parameter
dOPTIONS = odeset('RelTol', 1e-013, 'AbsTol', 1e-013);                      % Default integrator options
x0 = [1.00757; 0; 0.00257098670693846; 0; 0.0122713318389655; 0];           % Initial guess for Halo Orbit (hard-coded)
% x0 = richardson_plot(mu, 1, .28, 1);                                      % Initial guess for Halo Orbit (from Richardson approximation)
cnt = 0;

initConds = zeros(6, 1000);
perturbedConds = zeros(6, 360*1000);
unperturbedConds = zeros(6, 1000);
initTime = zeros(1, 1000);

%% Computation

x01 = [0.991650569223169;0;-0.000824723971902957;0;-0.00988275607142619;0]; % Already know at what point things get weird;
                                                                            % defined here. This should probably be changed

x0step = (x0(1)-x01(1))/1000;                                               % As above - defining that we want 1000 orbits between
                                                                            % the two extremes
x = x01;
                                                                           

for i = 1:1000                                                              % Existing functonality; continue orbit 1000 times
    
	[x, tf, success] = differentialCorrector(x0, mu);   
    % Differential correct the guess to ensure closure
    
    jacobiConst = jacobiConstant(x, 0., mu);
    
    if (jacobiConst <= 3.0007982727 && jacobiConst >= 3.0000030032)
        
        cnt = cnt + 1;
        
        initConds(:, cnt) = x;                                                    % Store the initial conditions and time for this orbit
        initTime(:, cnt) = tf;
        
        fprintf("Storing condition #%d\n", cnt);
        
%         [y2, y] = manifold(x,tf,mu);                                                 % Compute the conditions of the manifold at the |pi/8| section
%         perturbedConds(:,(cnt-1)*360+1:cnt*360) = y2(:,:); 
%         unperturbedConds(:,(cnt-1)*360+1:cnt*360) = y(:,1:6)';
                                                                                  % Continue the orbit
    else
        
        fprintf("J out of range. Skipping...\n");
        
    end
    
    initConds(:, i) = x;                                                    % Store the initial conditions and time for this orbit
    initTime(:, i) = tf;
    
%     y2 = manifold(x,tf,mu);                                                 % Compute the conditions of the manifold at the |pi/8| section
%     endConds(:,(i-1)*360+1:i*360) = y2(:,:)';                               % Store end conditions                                 
    x0 = x-[x0step;0;0;0;0;0];                                              % Continue the orbit
    fprintf('Computing orbit state for orbit number #%d\n', i);
    fprintf('%d\n', i);
    
    jacobiConst = jacobiConstant(x, 1, mu);
    fprintf("%12.10f", jacobiConst)
    if (jacobiConst >= 3.0008189806 || jacobiConst <= 3.0004448196)
        
        error("Jacobi constant out of range");
        
    end
   
end
    
fprintf('Done!');
save endConds2.mat endConds
save initTime2.mat initTime
save initConds2.mat initConds

%% Dependency subroutines

%
%	Generates EOM, using variational equations & STM
%

function xdot = cr3bp(t, x, mu)

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

%
%	Differential correction subroutine for periodicity checks
%

function [x, tf, success] = differentialCorrector(x0, mu)
	
    tolerance = 1e-12;                                                      % Set tolerance for periodicity checks					
	
    eventsOptions = odeset('Events', @findCross, ...
	'AbsTol', 1e-13, 'RelTol', 1e-13);                                      % Integrator options including events functions
	
    dimensions = 6;                                                         % Dimensions of phase space
	attemptNumber = 0;                                                      % Attempt number of differential correction sequence
	attemptNumberMax = 200;                                                 % Iteration controller
    success = false;                                                        % Iteration success flag
    
    while(attemptNumber <= attemptNumberMax)
        
        fprintf('Attempt number: %d', attemptNumber);
		
        phi0 = eye(6);                                                      % Initialise STM to identity matrix
		
        y0 = [x0; reshape(phi0, [], 1)];                                    % Assemble full initial state vector
		
        [t, y] = ode45(@cr3bp, ...
                      [0 10], y0, eventsOptions, mu);                       % Integrate to x-axis crossing
		
        tf = t(end,1);                                                      % Extract time of crossing
		xf = y(end, 1:dimensions);                                          % Extract s/c state at crossing
		
        if ((abs(xf(4)) < tolerance) && (abs(xf(6)) < tolerance))           % Periodicity check; if yes, break
			
            success = true;                                     
			break;
		
        end
        
		phi = reshape(y(end, dimensions+1:end), dimensions, []);            % Recover the state transition matrix
		
        % If we're here, then we were unsuccessful; begin computing
        % correctons
        
        % Compute accelerations...
        
       	r1 = sqrt((xf(1)+mu)^2 + xf(2)^2 + xf(3)^2);                        % Distance to primary
    	r2 = sqrt((xf(2)-1+mu)^2 + xf(2)^2+xf(3)^2);                        % Distance to secondary
		
        coeff1 = (1-mu)/r1^3;                                               % Define coefficients to make our life easier
		coeff2 = mu/r2^3;			 
        
		xddot = 2*xf(5)+xf(1)-...
                coeff1*(xf(1)+mu)-coeff2*(xf(1)-1+mu);                      % x-accel using potential function
		
        yddot = -2*xf(4)+xf(2)-coeff1*xf(2)...                  
                -coeff2*xf(2);                                              % y-accel; unneeded but kept in case
		
        zddot = -coeff1*xf(3)-coeff2*xf(3);                                 % z-accel

        fprintf('x error: %12.12f\t y error: %12.12f\n', ...
            xf(4)-tolerance, xf(6)-tolerance);
        
        % Construct the corrections matrix
        
        M1 = [phi(4,3) phi(4,5); ...                                        % z0 & ydot update (keep x0 fixed):
			 phi(6,3) phi(6,5)];

        M2 = -1/xf(5)*[xddot; zddot]*[phi(2,3) phi(2,5)];                   % Define second half separately; keep things clean
        
        M = M1 + M2;
			
        dzdydot0 = M\[-xf(4); -xf(6)];                                      % Solve with -xdot and -zdot because we only integrated to T/2:
                                                                            % \ => inverse(M) * ...
        x0 = x0 + [0; 0; dzdydot0(1); ...
                   0; dzdydot0(2); 0];                                      % Correct initial conditions and shoot again:
		    
        attemptNumber = attemptNumber + 1;                                  % Increment attempt number
	
    end
        
        x = x0;                                                             % Update output argument
        tf = 2.0*tf;

end

%
%	Events function - detects y-axis crossings
%

function [value, isterminal, direction] = findCross(t, x, mu)
	
    dimensions = 6;
	state = x(1:dimensions);
	
    value = state(2);
	isterminal = 1;
	direction = 0;

end
