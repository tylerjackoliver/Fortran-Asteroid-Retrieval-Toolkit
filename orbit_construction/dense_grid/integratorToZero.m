planeCondsZero = zeros(length(planeCondsDense)*1000, 6);

OPTIONS_events = odeset('RelTol', 1e-013, 'AbsTol', 1e-013, 'Events', 'on');

for i = 1:length(planeCondsDense)
    
    [t, x] = ode45('threeBodyEOM', [0 -100], planeCondsDense(i,:), OPTIONS_events);
    
    indices = linspace(1, length(x), 1000); % Generate 1000 points
    
    for j = 1:length(indices)
        
        temp = x(j, :);
        planeCondsZero((i-1)*1000+j,:) = x(floor(indices(j)), :);
        
    end
    
    fprintf("Completed state %d\n", i);
    
end

dlmwrite('2019-11-20_L2PlanarBackCondsSynodic.csv', planeCondsZero, 'delimiter', ',', 'precision', '%16.12f');