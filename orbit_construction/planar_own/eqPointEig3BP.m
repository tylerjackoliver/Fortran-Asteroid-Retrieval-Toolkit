function [Es,Eu,Ec,Vs,Vu,Vc] = eqPointEig3BP(ep,param)

%        [Es,Eu,Ec,Vs,Vu,Vc] = eqPointEig3BP(ep,mu);
%
% Find all the eigenvectors locally spanning the 3 subspaces for the phase
% space in an infinitesimal region around an equilibrium point 
%
% Our convention is to give the +y directed column eigenvectors
%
% NOTE: Uses Matlab's built-in numerical algorithms to solve the eigenvalue
%       equation which arises
%
% Shane Ross (revised 2.18.04)

Df = Dfmatrix3BP(ep,param);	 	% constant coefficient matrix at equil.

[Es,Eu,Ec,Vs,Vu,Vc] = eigGet(Df,0);	% find the eigenvectors
					% i.e., local manifold approximations 

% give the +y directed column vectors 
if Vs(2)<0, Vs=-Vs; end
if Vu(2)<0, Vu=-Vu; end
