% Creates Mach table
gamma = 5/3;
% First, P-M function
a = fopen('tabData.txt', 'w+');
fprintf(a, 'M,nu,mu \n');
for M = 1:0.01:3
    PMterm1 = sqrt((gamma-1)/(gamma+1)*(M^2-1));
    PMterm2 = sqrt(M^2-1);
    nuRad = sqrt((gamma+1)/(gamma-1))*atan(PMterm1)-atan(PMterm2);
    nuDeg = rad2deg(nuRad);
    waveAngleDeg = rad2deg(asin(1/M));
    fprintf(a, '%1.2f, %2.2f, %2.2f\n', M, nuDeg, waveAngleDeg);
end
fclose(a);