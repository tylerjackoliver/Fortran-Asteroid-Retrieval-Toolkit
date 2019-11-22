planeConds = zeros(length(planeCondsPruned*100), 6);

OPTIONS_events = odeset('RelTol', 1e-013, 'AbsTol', 1e-013, 'Events', 'on');

for i = 1:length(planeCondsPruned)
    
    [t, x] = ode45('threeBodyEOM', [0 -100], planeCondsPruned(:,i)', OPTIONS_events);
    
    planeConds(i, :) = x(end, 1:length100:end);
    
    fprintf("Completed state %d\n", i);
    
end