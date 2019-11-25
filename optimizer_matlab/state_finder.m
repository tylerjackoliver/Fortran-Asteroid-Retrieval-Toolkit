function [state_out, time_count, time_lower, time_upper] = ...
        state_finder(targ_can)
    
    %% Load the required SPICE kernels
    
    cspice_furnsh('de414.bsp');
    cspice_furnsh('naif0008.tls');
    
    %% Variable initialisation
    
    epoch = 0.0; % Start time
    epoch_counter = 0.0; % Epoch loop variable
    ang_can = 0.0; % Angle of candidate w.r.t. Earth
    ang_ear = 0.0; % Angle of the Earth w.r.t. Sun
    tau_can = 0.0; % Orbital period of the candidate
    a_can = 0.0;   % SMA of the candidate
    epoch_upper = 0.0; % Upper bound for the epoch
    last_checked = 0.0; % Time the state was last checked
    dum = 0.0;          % Dummy variable
    distance = 0.0;     % Distance between candidate and the Earth
    max_distance = 0.0; % Maximum distance between the candidate and the Earth
    state_can = zeros(6, 1); % State of the candidate
    state_ear = zeros(6, 1); % State of the Earth
    oe_can = zeros(6,1); % OE vector
    GMSun = 1.32712440018e11; %km^3/s^2
    iState = 0; % Program execution status flag
    
    %% SPICE variables
    
    abcorr = "NONE";
    obs = "Sun";
    targ_ear = "Earth";
    coord = "ECLIPJ2000";
    epoch_str = "Jan 1, 2020 00:00";
    epoch_upper_str = "Jan 1, 2050 00:00";
    
    %% Program execution
    
    % Load the corresponding SPICE kernel
    
    cspice_furnsh(sprintf('targ_can%s', '.bsp'));
    
    % Get the epoch in terms of ephemeris seconds (et)
    
    epoch = cspice_str2et(epoch_str);
    epoch_upper = cspice_str2et(epoch_upper_str);
    
    % Initialise last_checked to this time
    
    last_checked = epoch;
    
    % Loop through to determine coarse solution to problem

    for epoch_counter = epoch:86400.:epoch_upper
        
        state_can = cspice_spkezr(targ_can, epoch_counter, coord, ...
            abcorr, obs);
        state_ear = cspice_spkezr(targ_ear, epoch_counter, coord, ...
            abcorr, obs);
        
        % Get the orbital elements of the asteroid (to compute orbital
        % period)
        
        oe_can = cspice_oscelt(state_can, epoch_counter, GMSun);
        
        % Get angle between Earth and candidate
        
        ang_can = atan2(state_can(2), state_can(1));
        ang_ear = atan2(state_ear(2), state_ear(1));
        
        % Calculate SMA from rp, e
        
        a_can = oe_can(1) / (1.0 - oe_can(2));
        
        % Get period of candidate
        
        tau_can = 2 * pi * sqrt(a_can ^ 3 / GMSun);
        
        if (abs(ang_ear - ang_can) < pi / 8)
            
            last_checked = epoch_counter;
            
        end
        
        if ( (epoch_counter - last_checked) > tau_can)
            
            iState = 1;
            break;
            
        end
        
    end
    
    if (iState == 0)
        
        error('State finding solution cannot be found');
        
    end
    
    % If above was successful (iState = 1), then we can refine our solution
    
    time_lower = last_checked;
    time_upper = last_checked + tau_can;
    
    for epoch_counter = time_lower:16600:time_upper
        
        state_can = cspice_spkezr(targ_can, epoch_counter, coord, ...
            abcorr, obs);
        state_ear = cspice_spkezr(targ_ear, epoch_counter, coord, ...
            abcorr, obs);
        
        distance = norm( state_can - state_ear );
        
        if distance > max_distance
            
            max_distance = distance;
            state_out = state_can;
            time_count = epoch_counter;
            
        end
        
    end
    
    
    
end

