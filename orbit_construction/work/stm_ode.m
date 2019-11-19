function xdot = stm_ode(t, x, mu) 

%
%	Generates EOM, using variational equations & STM
%

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

%% Dependency subroutines

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