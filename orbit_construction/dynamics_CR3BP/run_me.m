% Script to compute and plot a family of periodic orbits and manifolds
%
% See chapter 4 of 
% Dynamical Systems, the Three-Body Problem and Space Mission Design 
% W.S. Koon, M.W. Lo, J.E. Marsden and S.D. Ross 
% ISBN 978-0-615-24095-4 
% Marsden Books, 2008.
%
% http://www.shaneross.com/books


% set parameters for algorithm

param = 3.00653E-06 ; % parameter of system (Sun-Earth)
%param = 1.E-4 ; % parameter of system (random)
% param = 1.2E-02 ; % parameter of system (close to Earth-Moon)
eqNum = 1;         % number of equilibrium point (L1 or L2?)
nFam  = 10;         % number of desired periodic orbits in family


gam = (param/3)^(1/3) ; % rough estimate of distance between secondary and eq pt


% first two amplitudes for continuation procedure to get p.o. family

Ax1  = 1e-3*gam    ; % initial amplitude (1 of 2)
Ax2  = 2*Ax1        ; % initial amplitude (2 of 2)


% get the initial conditions and periods for a family of periodic orbits

[x0po,T] = poFamGet3BP(param,eqNum,Ax1,Ax2,nFam) ;


% set accuracy for integration of individual trajectories on tube

%OPTIONS = odeset('RelTol',3e-6,'AbsTol',1e-6); % lowest accuracy
%OPTIONS = odeset('RelTol',3e-8,'AbsTol',1e-8); % med accuracy
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy


% plot the family of periodic orbits from the initial conditions

for k = nFam:-1:1,
  x0 = x0po(k,:) ;
  tb = 0 ;
  tf = T(k)/2 ;
  [x,t] = trajGet3BP(x0,tb,tf,param,OPTIONS) ;
  plot(x(:,1),x(:,2),'k',x(:,1),-x(:,2),'k');
  pause(0.01) ;
  hold on
  axis equal
end


% get and plot unstable manifold (negative branch) of a medium-sized p.o.
% specifically, plot 40 (roughly) equally spaced trajectories on this `tube' manifold
% whose initial conditions are displaced 1.e-6 from the p.o. in phase space

%%

% endConds_l1 = [];
% 
% for nMed = 1:length(x0po);
%     
%     [~, x0W] = poManiLocal3BP(x0po(nMed,:),T(nMed),0,-1,-1,1e-6,param,2*T(nMed),40);
%     
%     for j = 1:length(x0W)
%         
%         [t,x] = ode45('pcr3bp_2', [0 -10], x0W(j,:), odeset('AbsTol', 1e-012, 'RelTol', 1e-012, 'Events', 'on'));
%         endConds_l1 = [x(end,1) x(end,2) 0 x(end,3) x(end,4) 0; endConds_l1];
%     
%     end
%     
%     fprintf('Computing manifold conditions for orbit #%d\n', nMed);
%     
% end

% nMed = 4 ;
% [xW,x0W] = poManiLocal3BP(x0po(nMed,:),T(nMed),0,-1,-1,1e-6,param,2*T(nMed),40);

% [xW,x0W] =
% poManiLocal3BP(x0po(nMed,:),T(nMed),0,1,1,1e-6,param,2*T(nMed),40);