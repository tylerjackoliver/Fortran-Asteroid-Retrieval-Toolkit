jacobis = zeros(length(elts), 1);

for i = 1:length(jacobis)
    
    a = elts(1,i)/(1-elts(2,i)); % rp -> a
    jacobis(i) = 1/a + 2 * sqrt(a * (1-elts(2,i)^2)) * cos(elts(3,i));
    
end

fprintf("The maximum jacobi energy is: %16.12f\n", max(jacobis));
fprintf("The minimum jacobi energy is: %16.12f\n", min(jacobis));
fprintf("The minimum periapsis radius is: %16.12f\n", min(elts(1,:)));
fprintf("The maximum periapsis radius is: %16.12f\n", max(elts(1,:)));
