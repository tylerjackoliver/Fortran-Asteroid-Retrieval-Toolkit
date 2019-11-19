planeConds = zeros(length(perturbedConds), 6);

OPTIONS_events = odeset('RelTol', 1e-012, 'AbsTol', 1e-012, 'Events', 'on');

for i = 1:length(perturbedConds)
    
    [t, x] = ode45('threeBodyEOM', [0 -100], perturbedConds(:,i)', OPTIONS_events);
    
    planeConds(i, :) = x(end, :);
    
    fprintf("Completed state %d\n", i);
    
end