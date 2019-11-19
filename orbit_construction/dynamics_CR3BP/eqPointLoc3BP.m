function eqPos = eqPointLoc3BP(param,eqNum)

% 	 eqPos = eqPointLoc3BP(mu,eqNum)
%
% output:
% eqPos = position [x y] of equilibrium point on x-axis
%
% input:
% mu    = mass parameter of system
% eqNum = the number of the equilibrium point of interest
% 
% Shane Ross (revised 2.03.04)

mu = param ; % some parameter of system (mu)

mu1 = 1-mu ;
mu2 =   mu ;

% work out numbering convention here
 
lpts=zeros(5,2);
LPTS = [(mu1-gamma3BP(mu,1)) (mu1+gamma3BP(mu,2)) (-mu2-gamma3BP(mu,3))];
lpts(1,1)=LPTS(1);
lpts(2,1)=LPTS(2);
lpts(3,1)=LPTS(3);
lpts(4,1)= .5-mu;
lpts(5,1)= .5-mu;
lpts(4,2)= .5*sqrt(3);
lpts(5,2)=-.5*sqrt(3);

if 	eqNum == 1, eqPos = lpts(1,:) ; % L1
elseif 	eqNum == 2, eqPos = lpts(2,:) ; % L2
elseif 	eqNum == 3, eqPos = lpts(3,:) ; % L3
elseif 	eqNum == 4, eqPos = lpts(4,:) ; % L4
elseif 	eqNum == 5, eqPos = lpts(5,:) ; % L5
end
