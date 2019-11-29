[~, toInt] = size(initConds);
doptions = odeset('RelTol', 1e-012, 'AbsTol', 1e-012, 'Events', 'on');

newTime = zeros(toInt,1);

for i = 1:toInt
    
    if (initConds(1,i) == 0)
        
        continue;
        
    end
    
    [t,x] = ode45('threeBodyEOM', [0 10], initConds(:,i), doptions);
    
    newTime(i) = t(end) * 2;
    
    fprintf("Finished #%d\n; t is %f", i, t(end));
    
end

function [value, isterminal, direction] = findCross(~, x)
	
    dimensions = 6;
	state = x(1:dimensions);
	
    value = state(2);
	isterminal = 1;
	direction = -1;

end