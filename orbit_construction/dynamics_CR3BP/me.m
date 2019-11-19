% Set of scripts to single-shoot some orbits in about a lagrangian point
% L_i
%
% TODO: HOW THE FUCK DID I DO THE LINEARISATION?!

massParameter = 0.0121;                                                     % System mass parameter
lPoint = 1;                                                                 % Desired L-point
numberOfFamilies = 10;                                                      % Number of families to plot
lDist = (massParameter/3)^(1/3);                                            % Approximate distance to Li from secondary (Ross et. al, 2008)
amplitude1 = 1e-02*lDist;
amplitude2 = 2*amplitude1;                                                  % Get the first two sets of amplitudes for continuation procedure
[x0, T] = generateInitialGuess(massParameter, lPoint, amplitude1, ...       % Generate initial guess using (Richardson, 1980)
    amplitude2, numberOfFamilies);
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14);			    % Set default integrator options
% for k = 1:numberOfFamilies,
%   initialConds = x0(k,:);						    % Extract initial conditions from the IC array
%   tb = 0;								    % Set the integration time to begin at x0
%   tf = T(k)/2;								    % Integrate half of the orbit
%   [x,t] = integrator(initialConds, tb, tf, massParameter, OPTIONS);	
%   plot(x(:,1),x(:,2),'k',x(:,1),-x(:,2),'b');				    % Plot result
%   hold on
%   axis equal
%   [lLoc] = lPointFinder(massParameter, lPoint);				    % Append lagrangian point to plot -- find location
%   plot(lLoc(1), lLoc(2), '-ro');
% end
initialConds = x0(1,:);
% figure();
tf = T(1);
tb = 0;
[x, t] = integrator(initialConds, tb, tf, massParameter, OPTIONS);
plot(x(:,1),x(:,2),'k');
axis equal;
[lLoc] = lPointFinder(massParameter, lPoint);
plot(lLoc(1), lLoc(2), 'MarkerSize', 10);
function [x0, T] = generateInitialGuess(massParameter, lPoint, amplitude1,...
    amplitude2, numberOfFamilies)
    dt = -1e-02;                                                            % Guess at change in period between orbits
    dimensions = 4;                                                         % Number of dimensions we're dealing with (x(dot), y(dot))
    x0 = zeros(numberOfFamilies, dimensions);                               % Initialise matrices for states and periods
    T = zeros(numberOfFamilies, 1);
    
    x0Guess1 = initialGuessGenerator(massParameter, lPoint, amplitude1);             % Send these off for an initial guess(timate)
    x0Guess2 = initialGuessGenerator(massParameter, lPoint, amplitude2);
    
    % The following methodology is shamelessly copied from Ross, 2008: use
    % dx/dy from two different orbits in the family to seed the next set
       
    familyNumber = 1;
    disp(sprintf('Getting orbit family: %d', familyNumber));
    [x01, t1] = differentialCorrection(x0Guess1, massParameter);            % Send off the initial conditions to get refined guess

    familyNumber = 2;
    disp(sprintf('Getting orbit family: %d', familyNumber));    
    [x02, t2] = differentialCorrection(x0Guess2, massParameter);
    
    x0(1:2, 1:dimensions) = [x01(:)'; x02(:)'];                                      % Append the 'good' initial conditions to the matrix of init. conds.
    T(1:2) = [2*t1; 2*t2];                                                  % Transpose to get 4x1 -> 1x4
    
    for familyNumber = 3:numberOfFamilies                                   % Gotta catch 'em all!
         disp(sprintf('Getting orbit family: %d', familyNumber));
         dx = x0(familyNumber-1, 1) - x0(familyNumber-2, 1);                % Compute the deltas between the last two calc'd orbits
         dydot = x0(familyNumber-1, 4) - x0(familyNumber-2, 4);
         dt = T(familyNumber-1) - T(familyNumber-2);
         x0Guess = [(x0(familyNumber-1,1)+dx) 0 0 (x0(familyNumber-1,...
             4)+dydot)];
         [tempState, tempTime] = differentialCorrection(x0Guess, ...
             massParameter);
         x0(familyNumber, 1:dimensions) = tempState;
         T(familyNumber) = 2*tempTime;                                      % we only integrate 1/2, so x2!
    end
end


function [gamma] = pointDistance(lPoint, massParameter)                     % Find distance between lagrange and smaller primary
    m1 = massParameter;                                                     % Dependecy function for locating position of lPoint
    m2 = 1-massParameter;
    
    % Solve the quintic polynomial found in (Szebhely, 1967)
    l1Poly = [1 (m1-3) (3-2*m1) -m1 2*m1 -m1];
    l2Poly = [1 (3-m1) (3-2*m1) -m1 -2*m1 -m1];
    l3Poly = [1 (2+m1) (1+2*m1) -m2 -2*m2 -m2];
    
    l1Roots = roots(l1Poly);
    l2Roots = roots(l2Poly);
    l3Roots = roots(l3Poly);
    
    % Get rid of the complex pairs
    for k=1:5
        if isreal(l1Roots(k)) gammas(1)=l1Roots(k); end
        if isreal(l2Roots(k)) gammas(2)=l2Roots(k); end
        if isreal(l3Roots(k)) gammas(3)=l3Roots(k); end
    end
    gamma = gammas(lPoint);
end

function lPointPos = lPointFinder(massParameter, lPoint)
    m1 = 1-massParameter;
    m2 = massParameter;
    lPoints = zeros(5, 2);                                                  % Set up position matrix in planar case
    tempLPoints = [(m1 - pointDistance(1, massParameter)) ...               % Get the position of the three collienar points
        (m1 + pointDistance(2, massParameter)) (-m2...
        -pointDistance(3, massParameter))];
    lPoints(1, 1) = tempLPoints(1);                                         % Three collinear points are at y=0, so no need to assign from zeros matrix
    lPoints(2, 1) = tempLPoints(2);
    lPoints(3, 1) = tempLPoints(3);
    lPoints(4, 1) = 0.5-massParameter;                                      % X-coord of equilateral solutions
    lPoints(5, 1) = 0.5-massParameter;
    lPoints(4, 2) = 0.5*sqrt(3);                                            % Y-coord of equilateral solutions
    lPoints(5, 2) = -0.5*sqrt(3);
    
    lPointPos = lPoints(lPoint, :);                                         % Return desired position 
end

function jacobian = jacobianMatrix(lPos, massParameter)
    m1 = 1-massParameter;
    x = lPos;
    m2 = massParameter;
    r2= (x(1)+m2)^2+x(2)^2;      % r: distance to m1, LARGER MASS
    R2= (x(1)-m1)^2+x(2)^2;      % R: distance to m2, smaller mass
    % Calculate powers to ease up our life a little
    r3= r2^1.5;
    r5= r2^2.5;
    R3= R2^1.5;
    R5= R2^2.5;
    uxx = -1+(m1/r3)*(1-(3*(x(1)+m2)^2/r2))+(m2/R3)*(1-(3*(x(1)-m1)^2 ... % Compute the derivatives of the potential U
    /R2));
    uyy = -1+(m1/r3)*(1-(3*x(2)^2/r2))+(m2/R3)*(1-(3* x(2)^2/R2));
    uxy =   -(m1/r5)*3*x(2)*(x(1)+m2)-(m2/R5)*3*x(2)*(x(1)-m1);  
    jacobian = [0 0 1 0; 0 0 0 1; -uxx -uyy 0 2; -uxy -uyy -2 0];           % Arrange into restricted jacobian for planar case
end

function [eValStable, eValUnstable, eValCenter] = getEigenvalues(matrix)                                             % Grab eigenvalues & their stability from jacobian matrix around eq. points
    eValStable = [];
    eValUnstable = [];
    eValCenter = [];
    [height, width] = size(matrix);                                         % We need to analyse each bit of the returned matrix, so grab the dimensions
    [values, directions] = eig(matrix);                                     % Compute matrix eigenvals/vecs (not bothered about vec atm)
    stableCounter = 0;                                                      % We will categorise the nature of the eigenvalue, so initialise some counters
    unstableCounter = 0;
    centerCounter = 0;
    directions = cleanUpMatrix(directions);
    for i = 1:height
        if real(directions(i, i)) < 0                                       % Measure stability using sign
            stableCounter = stableCounter + 1;                              % Increment counter and add to array
            eValStable = directions(i, i);
        elseif real(directions(i, i)) > 0
            unstableCounter = unstableCounter + 1;
            eValUnstable = directions(i, i);
        else
            centerCounter = centerCounter + 1;
            eValCenter = directions(i, i);
        end
    end
end

function [eValStable, eValUnstable, eValCenter] = lPointEigenvalues(lPos,...
    massParameter)                                                           % Catch-all function to generate the eigvenVals about a point
    jacobian = jacobianMatrix(lPos, massParameter);
    [eValStable, eValUnstable, eValCenter] = getEigenvalues(jacobian);
end

function [x0Guess, tGuess] = initialGuessGenerator(massParameter,...
    lPoint, amplitude)                                                      % Generate an initial guess using (Richardson, 1980; Halo Orbit Formulation for the ISEE-3 mission)
    lPos = [lPointFinder(massParameter, lPoint) 0 0];                       % Set up position of collinear lagrange point
    [eValStable, eValUnstable, eValCenter] = lPointEigenvalues(lPos,...     % Get the eigenvalues for seeding initial guess from center eval
    massParameter);
    centerEigenvalueFrequency = abs(imag(eValCenter(1)));                   % Imaginary part of the center eigenvalue is the frequency
    gamma = pointDistance(lPoint, massParameter);                           % Get distance to L_i (Szebhely, 1967)
    l = centerEigenvalueFrequency(1);                                     % We want the frequency of the center eigenvalue, to avoid perturbing direction of propagation
    c2 = gamma^(-3)*(massParameter + (-1)^2*(1-massParameter)*...
        (gamma/(1-gamma))^(2+1));
    k = (l^2+1+2*c2)/(2*l);
    x0Guess = zeros(4, 1);                                                  % Initialise state matrix
    x0Guess(1) = lPos(1) - amplitude;                                       % Set up x guess
    x0Guess(4) = lPos(4) + amplitude*k*l;                                   % Set up y-dot guess
    tGuess = 2*pi/l;                                                        % Guess the time
end

function [x0, time1] = differentialCorrection(x0, massParameter)                % SINGLE SHOOT THE SHIT OUT OF THIS ORBIT
    MAXdxdot1 = 1e-010;
    RelTol = 3.e-14;
    AbsTol = 1.e-16;                                                         % Integration options
    maximumAttempts = 100;
    xdot1 = 1;
    attempt = 0;
    
    while abs(xdot1) > MAXdxdot1
        if attempt > maximumAttempts                                            % Computation handling
            disp('Oh shit, too many iterations!');
            break
        end
        model = 'cr3bp';                                                        % Point MATLAB to the governing equations
        time = [0 10];                                                          % Arbritrary time span for the integration
        options = odeset('RelTol', RelTol, 'AbsTol', AbsTol, 'Events', 'on');
        [temp1, temp2, time1, xx, temp3] = ode113(model, time, x0, options);      % Integrate until x-axis crossing
        x1 = xx(1);                                            % Extract state at x-axis crossing
        y1 = xx(2);
        xdot1 = xx(3);
        ydot1 = xx(4);
        if time1(1) < 1e-01
            time1 = time1(2);
        end
        options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
        [x, temp1, phi_t1, temp2] = stateTransitionMatrix(x0, time1, ...           % Propagate the STM over the orbit
            massParameter, options);
        attempt = attempt + 1;
        disp(fprintf('Differential correction: hoping this works number %d',...
            attempt));
	plot(x(:,1),x(:,2),'r-');
	hold on;
	m = length(x);
	plot(x(1,1),x(1,2));
	plot(x(m,1), x(1,2));
	pause(0.01);
   
        % Compute change in state for next attempt at correction
        m1 = 1-massParameter;                                                   % normalised masses of primaries in cr3bp
        m2 = massParameter;
        r3 = ((x1+m2)^2+y1^2)^1.5;                                              % r: distance to m1
        R3 = ((x1-m1)^2+y1^2)^1.5;                                              % R: distance to m2
        u_x = -x1+m1*(x1+m2)/r3 + m2*(x1-m1)/R3;                               % x-deriv. of the potential gunction
        xAccel = 2*ydot1-u_x;                                                   % Compute the change in xAccel
        yAccelDelta(attempt) = ...
            (1/(phi_t1(3,4)-phi_t1(2,4)*(1/ydot1)*xAccel))*xdot1;               % Compute the change in yAccel                               % Apply the changes
       if massParameter > 1e-03
           %dampFactor = 1-0.5^attempt;
           dampFactor = 1;
       else
           dampFactor = 1;
       end
        x0(4) = x0(4) - dampFactor*yAccelDelta(attempt);
    end
    disp('FUCK YEA, THAT WORKED');
end

function [x, time, phi, PHI] = stateTransitionMatrix(x0, time, ...
    massParameter, options, fixedStepFlag)
    dimension = 4;                                                          % Dimension of the problem; see above
    model = 'varEqs3BP';                                         % Point MATLAB to the variational equations model
    if nargin < 5                                                           % Argument control for fixedStepFlag
        fixed_step = 0 ;
        if nargin < 4
            OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14);                % integration options, if not specified
        end
    end
    
    if fixed_step == 0
        time1 = [0 time];
    else
        time1 = 0:time/(fixed_step-1):time;
    end
    PHI_0(1:dimension^2) = reshape(eye(dimension), dimension^2, 1);         % Initialise stm to (16x1) identity matrix
    PHI_0(1+dimension^2:dimension+dimension^2) = x0;                         % Initialise the end of PHI to be the state
    [time2, PHI] = ode113(model, time1, PHI_0, options);                       % Integrate the STM over the full orbit
    x = PHI(:,1+dimension^2:dimension+dimension^2); 		   % trajectory from time 0 to tf                       % Re-grab trajectory
    phi = reshape(PHI(length(time2), 1:dimension^2), dimension, dimension); % Grab the actual STM
end

 
function [x,t] = integrator(x0, tb, tf, massParameter, options)             % Time to integrate the FUCK out of the conditions, man
    model = 'cr3bp';
    options = odeset('RelTol',3e-10,'AbsTol',1e-12, 'Events', 'off');
    if nargin < 5                                                           % Option handling
        options = odeset('RelTol',3e-10,'AbsTol',1e-12, 'Events', 'off');
    end
    timetb = [-1e-012 tb];
    timetf = [0 tf];
    x = [];                                                                 % Initialise output array
    if tb == 0
        [t, x] = ode113(model, timetf, x0, options);
    elseif tf == 0
        [t, x] = ode113(model, timetb, x0, options);
    else
        [tt1, xx1] = ode113(model, timetb, x0, options);
        [tt2, xx2] = ode113(model, timetb, x0, options);
        x = [flipud(xx1);xx2];
        t = [flipud(tt1);tt2];
    end
end


function A = cleanUpMatrix(A)
    tolerance = 1e-014;
    for k = 1:length(A)
        for l = 1:length(A)
            if abs(real(A(k,1))) < tolerance
                A(k,l) = i*imag(A(k,l));
            end
            if abs(imag(A(k,l))) < tolerance
                A(k,l) = real(A(k,l));
            end
        end
    end
end

function manifoldIC = manifoldFinder(massParameter, x0, t0, dimensions, N)
    stm = eye(dimensions);                                                  % Initialise STM
    initialConds = [x0; stm];                                               % Initialise...initial conditions
    OPTIONS = odeset('RelTol', 1e-013, 'AbsTol', 1e-013);
    t0 = linspace(0, t0, N);
    [~, X] = ode45(@crtbp,  t0, initialConds, OPTIONS, massParameter);
    M = reshape(X(end, 5:end),6,[]);                                        % Extract monodromy matrix
    [monodromyEigVec, val] = eig(M);
    real = (imag(diag(val)) == 0);
    [~, i1] = max(abs(val));
    [~, i2] = min(abs(val));
    a = 1:6;
    stable = monodromyEigVec(:, a(i1));
    unstable = monodromyEigVec(:, a(i2));
    xs_in = zeros(dimensions, N);
    xs_out = xs_in;
    xu_in = xs_in;
    xu_out = xs_in;
    for i = 1:N
        Phi = reshape(X(i,5:end), []);
        x = X(i, 1:4);
        S = Phi*stable;
        U = Phi*unstable;
        xs_out(:, i) = x(:, i) + S.*10^(-4); 
        xs_in(:,i) = x(:, i) - S.*10^(-4);
        xu_out(:,i) = x(:,i) + U.*10^(-4);
        xu_in(:,i) = x(:,i) - U.*10^(-4);
    end
    manifoldIC = [xs_out xs_in xu_out xu_in];
end
    
end