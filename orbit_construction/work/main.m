
% Attempts a revised form of the differential correction routine for
% vertical lyapunov orbits

% Initialise variables

tolerance = 1e-010;
mu = 3.003458e-06;
maxAttempts = 200;

% Initialise structures

OPTIONS_events = odeset('RelTol', tolerance, 'AbsTol', tolerance, ...
                        'Events', 'off');
OPTIONS = odeset(OPTIONS_events, 'Events', 'off');

% Initial guess goes here

x0 = [1.00049 0 0.01897 0 -0.00740 0]';

x0step = 1e-05;                                                             % As above - defining that we want 1000 orbits between
                                                                            % the two extremes
endConds = zeros(6, 1000*360);                                              % Initialise output arrays (360 points per orbit * 1000 orbits)
initConds = zeros(6, 1000);
initTime = zeros(6, 1000);

for i = 1:1000                                                              % Existing functonality; continue orbit 1000 times

	[x,tf, success] = differentialCorrector(x0, mu);                        % Differential correct the guess to ensure closure
    
    if (not(success))
        
        error("Differential correction failed.");
        
    end
    
    initConds(:, i) = x;                                                    % Store the initial conditions and time for this orbit
    initTime(:, i) = tf;
    
    y2 = manifold(x,tf,mu);                                                 % Compute the conditions of the manifold at the |pi/8| section
    endConds(:,(i-1)*360+1:i*360) = y2(:,:)';                               % Store end conditions                                 
    x0 = x-[x0step;0;0;0;0;0];                                              % Continue the orbit
    fprintf('Computing orbit state for orbit number #%d\n', i);
    fprintf('%d\n', i);
    jacobiConst = jacobiConstant(x, mu);
    
    if (jacobiConst >= 3.0007982727 || jacobiConst <= 3.0000030032)
       
        error("Jacobi constant out of range");
        
    end

    hold on;
    
    [t, x] = ode45('threeBodyEOM', [0 tf], x, OPTIONS);
    plot3(x(:,1), x(:,2), x(:,3));
   
end
    
fprintf('Done!');
save endConds2.mat endConds
save initTime2.mat initTime
save initConds2.mat initConds

