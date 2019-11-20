[len,~] = size(endConds);

elements = [];

states = [];

mu = 1.0;

cspice_furnsh('naif0008.tls');

jd = cspice_str2et('01 Jan 2000 12:00');

cnt = 1;

times = [];

for i = jd:86400:(jd + 365.25 * 86400)
    
    a = cspice_et2utc(i, 'J', 8);
    mjd = str2double(a(4:end)) - 2400000.5;
    state = invRotationMatrix_final(endConds(55,:), mjd);                   % Rotate in at J2000
    states(cnt, :) = state;
    times = [times; mjd];
    elements(cnt, :) = cspice_oscelt(state', i, mu)'; 
    fprintf("Completed state #%d\n", cnt);
    cnt = cnt+1;
    
end

cspice_kclear;