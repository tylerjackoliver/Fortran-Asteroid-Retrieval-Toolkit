endConds_l1 = [];

for nMed = 1:length(x0po);
    
    [~, x0W] = poManiLocal3BP(x0po(nMed,:),T(nMed),0,-1,-1,1e-6,param,2*T(nMed),40);
    
    for j = 1:length(x0W)
        
        [t,x] = ode45('pcr3bp_2', [0 -10], x0W(j,:), odeset('AbsTol', 1e-012, 'RelTol', 1e-012, 'Events', 'on'));
        endConds_l1 = [x(end,1) x(end,2) 0 x(end,3) x(end,4) 0; endConds_l1];
    
    end
    
    fprintf('Computing manifold conditions for orbit #%d\n', nMed);
    
end