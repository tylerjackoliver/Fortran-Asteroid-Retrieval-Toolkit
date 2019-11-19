for i = 1:359999
        
    if (abs(atan2(planeConds(i, 2), planeConds(i, 1)) - pi/8) > 1e-06)
        
        fprintf("First case! Difference is %f\n", abs(atan2(planeConds(i, 2), planeConds(i, 1)) - pi/8));
        
    end
       
end