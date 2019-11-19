function [x0po,T] = poFamGet3BP(mu,eqNum,Ax1,Ax2,nFam) ;

%        [x0po,T] = poFamGet3BP(mu,eqNum,Ax1,Ax2,nFam) ;
%
% Generate a family of periodic orbits (po) given a pair of 
% seed initial conditions and periods
%
%----------------------------------------------------------------------------
% CR3BP with the LARGER MASS, m1 to the left of the origin at (-mu,0)
% and m2, or the planet (ie. Earth), at (1 - mu, 0)
%
%                L4
%
% -L3----m1--+-------L1--m2--L2-
%
%                L5
%
% Shane Ross (revised 2.19.04)

param = mu ;

% delt = guessed change in period between successive orbits in family

delt = -1.e-2 ;   % <==== may need to be changed

N = 4 ; % dimension of phase space

x0po = zeros(nFam,N) ;
T    = zeros(nFam,1) ;

[x0poGuess1,TGuess1] = poGuessLinear3BP(param,eqNum,Ax1) ;
[x0poGuess2,TGuess2] = poGuessLinear3BP(param,eqNum,Ax2) ;

% Get the first two periodic orbit initial conditions

  iFam = 1 ;
  FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
  disp(FAMNUM) ;
  [x01,t1] = poDifCor3BP(x0poGuess1,param) ;                  

  iFam = 2 ;
  FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
  disp(FAMNUM) ;
  [x02,t2] = poDifCor3BP(x0poGuess2,param) ;                  

x0po(1:2,1:N) = [x01(:)' ; x02(:)'] ;
T   (1:2)     = [2*t1    ; 2*t2   ] ;

for iFam = 3:nFam,

  FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
  disp(FAMNUM) ;

  dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
  dyp = x0po(iFam-1,4) - x0po(iFam-2,4) ;

  dt  = T(iFam-1) - T(iFam-2) ;

  x0po_g = [ (x0po(iFam-1,1) + dx) 0 0 (x0po(iFam-1,4) + dyp) ] ;
  t1po_g =   (T(iFam-1) + dt)/2 + delt ;
  
  % differential correction takes place in the following function
  
  [x0po_iFam,t1_iFam] = poDifCor3BP(x0po_g,mu) ; 

  x0po(iFam,1:N) = x0po_iFam ;
  T   (iFam)     = 2*t1_iFam ;	

  if mod(iFam,10) == 0,
    dum = [x0po T] ;
    save x0po_T.dat -ascii -double dum
  end

end

dum = [x0po T] ;
save x0po_T.dat -ascii -double dum
