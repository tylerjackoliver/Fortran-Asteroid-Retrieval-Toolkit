function [x,t] = trajGet3BP_2(x0,tb,tf,mu,OPTIONS);

%        [x,t] = trajGet3BP(x0,tb,tf,mu,OPTIONS);
%
% This is an integrator for a point in phase space
%
%        d x(t)
%        ------ =  f(x)
%          dt
%
%-----------------------------------------------------------
%
% Shane Ross (revised 10.13.03)
global param

param = mu ;
MODEL = 'pcr3bp_2';

if nargin==4;
   OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12, 'Events', 'on');  % lower accuracy
%  OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-16);  % high accuracy
end

TSPANtb = [-1.e-12 tb ] ;
TSPANtf = [ 0      tf ] ;

x=[];

if tb==0,
  [t,x]     = ode113(MODEL,TSPANtf,x0,OPTIONS);
elseif tf==0,
  [t,x]     = ode113(MODEL,TSPANtb,x0,OPTIONS);
else,
  [tt1,xx1] = ode113(MODEL,TSPANtb,x0,OPTIONS);
  [tt2,xx2] = ode113(MODEL,TSPANtf,x0,OPTIONS);

  x=[flipud(xx1);xx2];
  t=[flipud(tt1);tt2];
end
