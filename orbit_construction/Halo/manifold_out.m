% File: manifold.m
% Author: Nick Truesdale
% Date: 12/12/2012
% 
% Description: This script creates invariant manifolds for CRTBP periodic
% orbits. 

% Housekeeping
function Y2_end = manifold(y0,t0)

% Constants
AU = 149597870.691; % km
year = 0.99997862*365.25*24*3600;
N = 360;  % number of trajectories to consider
dx = 500; % km, perturbation from Halo orbit
tshort = 1.35*pi;
tlong = 2*pi;

% CRTBP parameters
M_earth = 1;
M_sun = 328902.82113001; % Relative to earth moon bary
M_moon =  0.0123000569113856;

mu = M_earth/(M_sun + M_earth);
d = AU;
xpert = dx/d;

% Lagrange points
L = lagrange_points(mu);

% Initial conditions
orbit = 2;
% load init
% y0 = init(:,orbit);
% t0 = t0(orbit);
% y0 = [1.00765; 0; 0.00257098670693846; 0; 0.0122713318389655; 0];
% t0 = 3.09069728945353;

% y0 = [1.007650000000000 0 0.468449229815266 0 -0.113349292825486 0]';
% t0 = 6.281929581429136;

% y0 = [1.007580111481132; 0; -0.0027; 0; 0.012542853074521; 0];
% t0 = 3.08982093888028;

% y0 = [0.991975569223169;0;-0.001886564439553;0;-0.010972997563405;0];
% t0 = 3.0557;

% y0 = [0.996155555377277,0,-0.00504623798961148,0,-0.0232146192170136,0]';
% t0 = 2.933695679823386;

% y0 = [0.991975569223169;0;-0.001886564439553;0;-0.010972997563405;0];
% y0 = [0.991975569223169;0;0.001886564439553;0;-0.010972997563405;0];
% t0 = 1.527856078270380;

% y0 = [0.991975555377273,0,-0.00191718187217543,0,-0.0110295021073721,0]';
% t0 = 3.05553470727115;

% y0 = [0.991975569223169;0;-0.001886564439553;0;-0.010972997563405;0]
% 
% y0 = [1.007620000000000;0;0.002617959033317;0;0.012364637741252;0];
% t0 = 3.090168748897602;
% 
% y0 = [1.007570000000000;0;0.002693911126194;0;0.012519612580424;0];
% t0 = 3.089276602226325;
% 
% y0 = [1.007540000000000;0;0.002738163450494;0;0.012612282813562;0];
% t0 = 3.088734494959405;

% Get manifold initial conditions
[x1, x2, x3, x4, Y] = manifold_mono(mu, y0, t0, N, xpert);
xnom = Y(:,1:6)'; 

% Axes setup
d = 1 - L(1); % distance unit
ax1 = [1-1*d, 1+6*d, -0.5*d, 7*d, -0.01, 0.03];
ax2 = [1-0.5*d, 1+1.5*d, -d, d];
ax3 = [1-1*d, 1+6*d, -7*d, 0.5*d, -0.01, 0.03];
ax4 = [1-0.75*d, 1+1.25*d, -d, d];

%% Stable Manifold Plots
% Set up long figure
% scr = get(0,'ScreenSize');
% h1 = figure('OuterPosition', [scr(3)*0.5, 200, scr(4)*0.8, scr(4)*0.8]);
% h1 = figure;
hold on
% plot3(-mu, 0, 0, 'ok', 'markerfacecolor', 'y', 'markersize', 30)
% plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 22)
% grid on
% axis equal
% axis(ax1)
% view(-25, 16)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Stable Manifold')

OPTIONS_events = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', 'on');
Y2_end = [];
% Loop through perturbed initial conditions
for ii = 1:N    
    % Propogate orbit for one revolution (pi) and gather data at N points
%     opts = odeset('reltol', 1e-13, 'abstol', 1e-13);
    [~, Y1] = ode45('threeBodyEOM', [0 -3*tlong], x2(:,ii), OPTIONS_events);
%     [~, Y2,~,~,~] = ode45('threeBodyEOM', [0 -3*tlong], x1(:,ii), OPTIONS_events);
    Y2_end = [Y2_end; Y1(end,:)];
    
    % Long Plot
%     plot3(Y1(:,1), Y1(:,2), Y1(:,3), 'color', 'blue')
%     plot3(Y2(:,1), Y2(:,2), Y2(:,3), 'color', 'black')
    
end

% plot3(L(:,1), L(:,2), zeros(5,1), 'ok', 'markerfacecolor', 'r', 'markersize', 6)
%saveas(h1, '../writeup/images/stable0', 'png')
% axis(ax2)
view(0, 90)
%saveas(h1, '../writeup/images/stable1', 'png')

end

% %% Unstable Manifold
% % Set up figure
% % h2 = figure('OuterPosition', [scr(3)*0.5, 100, scr(4)*0.8, scr(4)*0.8]);
% h2 = figure;
% hold on
% plot3(-mu, 0, 0, 'ok', 'markerfacecolor', 'y', 'markersize', 30)
% plot3(1-mu, 0, 0, 'ok', 'markerfacecolor', 'b', 'markersize', 15)
% grid on
% axis equal
% axis(ax4)
% view(0, 90)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Unstable Manifold')
% view(205, 16)
% 
% % Nominal orbit
% t = linspace(0, tshort, 2000);
% a = 1:length(t);
% opts = odeset('reltol', 1e-13, 'abstol', 1e-13);
% [~, Y] = ode45(@three_body_ode, t, y0, opts, mu);
% 
% % Analysis Parameters
% mu_earth = 39860.4415;
% r_em = 384403; v_em = sqrt(mu_earth/r_em);
% r_geo = 35786; v_geo = sqrt(mu_earth/r_geo);
% earth = [1-mu; 0; 0];
% 
% dist = 15000/AU;
% geo = r_geo/AU;
% em = r_em/AU;  
% 
% % Preallocate
% nn = 1;
% mm = 1;
% 
% % Loop through perturbed initial conditions
% for ii = 1:N    
%     % Propagate orbit for one revolution (pi) and gather data at N points
%     [~, Y] = ode45(@three_body_ode, t, xnom(:,ii), opts, mu);
%     [~, Y1] = ode45(@three_body_ode, t, x3(:,ii), opts, mu);
%     [~, Y2] = ode45(@three_body_ode, [0 tlong], x4(:,ii), opts, mu);
%     
%     % Find time when satellite has "left" the orbit
%     dx = colnorm(Y1(:, 1:3)' - Y(:, 1:3)');
%     
%     ind1 = dx < dist;
%     t_left = t(ind1); t_left = t_left(end);
%     
%     % Find trajectories that approach Earth
%     dx = colnorm(Y1(:, 1:3)' - repmat(earth, 1, length(Y1))); 
%     ind2 = dx < em;
%     ind3 = dx < geo;
%     
%     if(sum(ind3))   
%         % Find closest approach
%         y = Y1(ind3, :); y = y(1,:);
%         tt = t(ind3);
%         t3(nn) = (tt(1) - t_left)/2/pi;
%         theta3(nn) = (ii-1)/(N-1)*360;
%         
%         % Find velocity at closest approach and delta V required to insert
%         v3(nn) = norm(y(4:6))*AU/year*2*pi;
%         dv3(nn) = sqrt((v3(nn)^2 - v_geo^2) * 2*mu_earth/r_geo);
%         
%         % Coloring and increment
%         C = 'k';
%         l = 2;
%         nn = nn + 1;
%     end
%         
%     if(sum(ind2))
%         % Find closest approach
%         y = Y1(ind2, :); y = y(1,:);
%         tt = t(ind2);
%         t2(mm) = (tt(1) - t_left)/2/pi;
%         
%         % Find velocity at closest approach and delta V required to insert
%         v2(mm) = norm(y(4:6))*AU/year*2*pi;
%         dv2(mm) = sqrt((v2(mm)^2 - v_em^2) * 2*mu_earth/r_em);
%         theta2(mm) = (ii-1)/(N-1)*360;
%         
%         % Coloring and increment
%         if(sum(ind3) == 0)
%             C = 'b'; 
%             l = 1;
%         end
%         mm = mm + 1;
%         
%     else
%         C = 'g';
%         l = 0.5;
%     end
%     
%     % Outward
%     plot3(Y2(:,1), Y2(:,2), Y2(:,3), 'color', [0.8, 0.95, 0.85])
%     
%     % Plot Inward
%     plot3(Y1(:,1), Y1(:,2), Y1(:,3), 'color', C, 'linewidth', l)
% end
% 
% % Complete figure
% plot3(L(:,1), L(:,2), zeros(5,1), 'ok', 'markerfacecolor', 'r', 'markersize', 6)
% % saveas(h2, '../writeup/images/unstable0', 'png')
% axis(ax4)
% view(0, 90)
% % saveas(h2, '../writeup/images/unstable1', 'png')
% 
% clearvars -except Y v2 v3 t2 t3 theta2 theta3
% 
% %% Plots for intersecting orbits
% h5 = figure;
% hold on
% plot(theta2, t2*365.25, 'b.', theta3, t3*365.25, 'k.')
% plot(theta2, t2*365.25, 'b:', theta3, t3*365.25, 'k:', 'linewidth', 2)
% xlabel('Angle around Halo Orbit (deg)')
% ylabel('Time of Flight (days)')
% title('Flight Time for Moon and Geosynchronous Transfers')
% legend('To Lunar Orbit', 'To Geosynch.', 'location', 'northeast')
% grid on
% 
% h6 = figure;
% hold on
% plot(theta2, v2, 'b.', theta3, v3, 'k.')
% plot(theta2, v2, 'b:', theta3, v3, 'k:', 'linewidth', 2)
% xlabel('Angle around Halo Orbit (deg)')
% ylabel('Insertion Delta V (km/s)')
% title('Delta V for Moon and Geosynchronous Transfers')
% legend('To Lunar Orbit', 'To Geosynch.', 'location', 'northeast')
% grid on
% 
% saveas(h5, '../writeup/images/time', 'png')
% saveas(h6, '../writeup/images/dv', 'png')