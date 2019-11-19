%
% Computes eigenvalues of the Monodromy matrix for some orbit y0 with
% period t0, and integrates the resulting perturbed states to the |pi/8|
% plane
%

function Y_end = manifold(y0,t0, mu)
    %% Initialisation

    N = 60;                                                                    % Number of discretisations to compute

    tshort = 1.35*pi;                                                           % Plot integration time; short
    tlong = 2*pi;                                                               % Plot integration time; long

    xpert = 1e-06;                                                              % Magnitude of the x-perturbation
    
    % Perturb manifold initial conditions

    [XS_in, XS_out, ~, ~, ~] = manifold_mono(mu, y0, t0, N, xpert);             % XS_in, XS_out, XU_in, XU_out, Monodromy matrix

    %% Computation

    OPTIONS_events = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, ...
                            'Events', 'on');                                    % Events -> stop at |pi/8| plane

    Y_end = zeros(N, 6);                                                       % Final conditions array

    hold on;                                                                    % If we want to plot

    %% Check eigenvalue direction
    
    % Modification 7/3/18 - for L2-based orbits, we want to perturb our
    % orbit in the direction s.t. the pertubation leads to a
    
    % Loop through perturbed initial conditions
    % We need to choose the correct manifold direction (want stable in
    % backwards time) - whether we perturb in or out. Thus, adjust as 
    % necessary below to get the correct direction (i.e. choose x1 or x2)
    % Therefore, compute a reduced set of the integration and check all is
    % in the right direction
    
    test_vals_1 = [];
    test_vals_2 = [];
    
    for ii = 1:N/12:N
        
        [~, Y1] = ode45('threeBodyEOM', [0 -3*pi], XS_in(:,ii), ...
                 OPTIONS_events);
             
        [~, Y2] = ode45('threeBodyEOM', [0 -3*pi], XS_out(:,ii), ...
                      OPTIONS_events);                                          % Compute points at intersection of |pi/8| plane
        
        test_vals_1(end+1) = Y1(end, 2);
        test_vals_2(end+1) = Y2(end, 2);
                  
    end
    
    % For L1 orbits, we desire all the points to be negative...
    
    if (test_vals_1(:) < 0)
        
        dset = XS_in;
        
    elseif (test_vals_2(:) < 0)
        
        dset = XS_out;
        
    else
        
        error("Neither directions for eigenvalue perturbations are OK!");
        
    end
    
    % Now integrate all of the 'correct' directions
    
    for ii = 1:N

        % Propagate orbit for one revolution (pi) and gather data at N 
        % points

        [~, Y1] = ode45('threeBodyEOM', [0 -10*tlong], dset(:,ii), ...
                         OPTIONS_events);

        Y_end(ii,:) = Y1(end,:);

        % For debugging - plot each of the directions and see which one is
        % correct

        plot3(Y1(:,1), Y1(:,2), Y1(:,3), 'color', 'blue');

    end

end

