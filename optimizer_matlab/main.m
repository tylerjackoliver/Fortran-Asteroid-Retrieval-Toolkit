% Set up and run the optimisation problem for asteroids and shit

cspice_furnsh('naif0008.tls');
transfer_epoch = cspice_str2et('23 Sep 2036 00:00');

% Set bounds

XL = [transfer_epoch * .99, 1507* 86400 * 0.8];
XU = [transfer_epoch * 1.01, 1507 * 86400 * 1.2];

% Define initial guess

X0 = [transfer_epoch, 1507 * 86400];

% Call fmincon, give it a whack

X = fmincon(@cost_function, X0, [], [], [], [], XL, XU);

% Call cost function one more time

cost_function(X)