function  gamma = gamma3BP(mu, Lpt)

% gamma = gamma3BP(mu, Lpt);
%
% Calculate ratio of libration point distance from closest primary to distance
% between two primaries, e.g. gamma(mu,2) = distance between M2 and L2
%
%----------------------------------------------------------------------
% CR3BP (Circular Restricted Three-Body [Gravitational] Problem)
% with the LARGER MASS, M1 (i.e., Sun) to the left of the origin at (-mu,0)
% and the smaller mass, M2, or the planet (i.e., Earth), is at (1-mu, 0)
%
%      (rotating frame)
%
%                 L4
%
% -L3------M1--+-----L1--M2--L2-
%
%                 L5
%			 |<-->|
%			  gamma
%
% Shane Ross (revised 1.29.04)

mu2 = 1 - mu;

%        x^5    x^4       x^3      x^2   x    const
poly1 = [1   -1*(3-mu)  (3-2*mu)  -mu   2*mu  -mu ];
poly2 = [1      (3-mu)  (3-2*mu)  -mu  -2*mu  -mu ];
poly3 = [1      (2+mu)  (1+2*mu)  -mu2 -2*mu2 -mu2];

% solve for roots of quintic polynomial

rt1 = roots(poly1);
rt2 = roots(poly2);
rt3 = roots(poly3);

% keep only the real root (there are also two complex pairs)

for k=1:5
        if isreal(rt1(k)) GAMMAS(1)=rt1(k); end
        if isreal(rt2(k)) GAMMAS(2)=rt2(k); end
        if isreal(rt3(k)) GAMMAS(3)=rt3(k); end
end

gamma = GAMMAS(Lpt) ;
