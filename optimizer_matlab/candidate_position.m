function state_out = candidate_position(t_epoch, state_epoch, state_in)

    % Since MATLAB does not have PROP2B, we need to rotate into
    % the three-body problem and then integrate forward
    
    au = 1.495978707e8;
    
    % Non-dimensionalise
    
    state_in(1:3) = state_in(1:3)/au;
    state_in(4:6) = state_in(4:6) * 365.25 * 86400 / (2 * pi * au);
    
    % Rotate
    
    state_cr3bp = globalTocr3bp(state_in, state_epoch);
    
    % Determine time delta in non-dim units
    
    delta_time = (t_epoch - state_epoch) * 2 * pi / (365.25 * 86400);
    
    % Integrate
    
    [~, x] = ode45('threeBodyEOM', [0 delta_time], state_cr3bp, ...
        odeset('RelTol', 1e-012, 'AbsTol', 1e-012));
    
    % Extract final state
    
    state_out_cr3bp = x(end, :);
    
    state_out = cr3bpToGlobal(state_out_cr3bp, t_epoch);
    
    % Dimensionalise
    
    state_out(1:3) = state_in(1:3)*au;
    state_out(4:6) = (state_in(4:6) / (365.25 * 86400)) * (2 * pi * au);

end