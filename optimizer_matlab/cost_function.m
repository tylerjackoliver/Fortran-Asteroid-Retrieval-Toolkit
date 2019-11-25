function min_vel = cost_function(x)
%COST_FUNCTION Cost function for asteroid retrieval

    %% Variable initialisation

    state_can = zeros(6,1);         % Asteroid candidate state
    state_can_orig = zeros(6,1);    % Un-rotated candidate state
    state_targ = zeros(6,1);        % Original state of the target
    state_rot  = zeros(6,1);        % Rotated state of the target
    transfer_epoch = 0.0;           % Epoch of the transfers
    transfer_time = 0.0;            % Transfer time
    state_epoch = 0.0;              % Epoch of the state
    time_lower = 0.0;               % Lower bound for transfer_epoch
    time_upper = 0.0;               % Upper bound for transfer_epoch
    tt = 0.0;                       % Transfer time
    vx1 = 0.0; vx2 = 0.0;           % x-velocities
    vy1 = 0.0; vy2 = 0.0;           % y-velocities
    vz1 = 0.0; vz2 = 0.0;           % z-velocities
    vtx = 0.0; vty = 0.0;           % Transfer velocities
    vtz = 0.0; vcx = 0.0;           % Transfer velocities
    vcy = 0.0; vcz = 0.0;           % Transfer velocities
    
    long_way = false;               % Boolean for Lambert solver
    run_ok   = false;               % Boolean for Lambert solver
    
    v1 = [];                        % Lambert arc velocities
    v2 = [];                        % Lambert arc velocities
    
    multi_rev = 4;                  % Number of revolutions on the Lambert
    
    targ_can = '3435539';           % Candidate to optimise
    
    %% Get the state of the candidate, and the epoch bounds
    
    [state_can_orig, state_epoch, ~, ~] = ...
        state_finder(targ_can);
    
    %% Extract problem parameters
    
    transfer_epoch = x(1);
    tt = x(2);
    
    %% Initialise velocities
    
    min_vel = 1.e6;
    
    %% Open the input file
    
    top_transfers = dlmread('2019-11-23_topTransfers50000.csv');
    
    %% Rotate candidate state from the state_epoch to the current time
    
    state_can = candidate_position(transfer_epoch, state_epoch, ...
        state_can_orig, state_can);
    
    %% Perform main optimisation loop
    
    for i = 1:length(top_transfers)
        
        state_targ = top_transfers(i, :);
        
        % The states above are defined at J2000. Therefore, apply the
        % relations from Sanchez et. al. to "propagate" the states forward
        
        state_rot = rotator(state_targ, transfer_epoch, tt);
        
        % Define the initial and final velocities of the target and the
        % candidate
        
        vcx = state_can(4); vcy = state_can(5); vcz = state_can(6);
        vtx = state_rot(4); vty = state_rot(5); vtz = state_rot(6);
        
        % Call the Lambert solver on state_rot (the target rotated to the
        % "correct" epoch). !!tt is in DAYS!!
        
        [v1, v2, ~, ~] = lambert(r1, r2, tf, multi_rev, GM_central);
        
        % Now iterate over v1, v2 to find transfer velocities
        
        for j = 1:(numel(v2)/3)
            
            vx1 = v1(j, 1);                                                                                                     % x-velocity at start of Lambert arc
            vy1 = v1(j, 2);                                                                                                     % y-velocity at start of Lambert arc
            vz1 = v1(j, 3);                                                                                                     % z-velocity at start of Lambert arc

            vx2 = v2(j, 1);                                                                                                     % x-velocity at end of Lambert arc
            vy2 = v2(j, 2);                                                                                                     % y-velocity at end of Lambert arc
            vz2 = v2(j, 3);                                                                                                     % z-velocity at end of Lambert arc

            transfer_v1  = ((vx1-vcx)^2 + (vy1-vcy)^2. + (vz1-vcz) ^ 2.)^.5;                                                    % deltaV from the candidate and the beginning of the lambert arc
            transfer_v2  = ((vx2-vtx)^2 + (vy2-vty)^2. + (vz2-vtz) ^ 2.)^.5;                                                    % deltaV from the target and the end of the lambert arc
            transfer_vel = transfer_v1 + transfer_v2;                                                                           % Total delta v is the sum of both; both > 0

            if transfer_vel < min_vel
                
                min_vel = transfer_vel;
                
            end
            
        end
            
    end
    
end

