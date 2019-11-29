planeCondsDense = zeros(length(toInt*100), 6);

OPTIONS_events = odeset('RelTol', 1e-013, 'AbsTol', 1e-013, 'Events', 'on');

for i = 1:length(toInt)
    
    [t, x] = ode45('threeBodyEOM', [0 -100], toInt(:,i)', OPTIONS_events);
    
    planeCondsDense(i, :) = x(end, :);
    
    fprintf("Completed state %d\n", i);
    
end