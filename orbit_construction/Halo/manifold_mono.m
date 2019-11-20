function [XS_in, XS_out, Y] = ...
          manifold_mono(mu, y0, t0, N, xpert)

%
% This function propagates the initial condition Y0 for one 
% period (up to T0) in the CRTBP frame defined by MU. The monodromy matrix 
% is determined, and its eigenvectors used to determine perturbations of 
% magnitude xpert at N points along the orbit.
%
% Each of the outputs are 6 by N arrays whose columns are state vectors.
% The XS terms are for the stable manifold and the XU terms are unstable. 
% The "in" and "out" arrays correspond to positive and negative
% perturbations.
%

    % Initial conditions - add STM

    stm = reshape(eye(6), 36, 1);
    y0 = [y0; stm];

    % Set up ODE45
    t = linspace(0, t0, N);                                                     % Gather data at N points
    opts = odeset('reltol', 1e-13, 'abstol', 1e-13);

    % Propagate orbit for one revolution
    [~, Y] = ode45(@stm_ode, t, y0, opts);                                  % Integrate with a special form of the EOM:
                                                                                % includes propagation of the STM

    % Monodromy matrix

    M = reshape(Y(end, 7:42), 6, 6);

    % Eigenvalue and eigenvector analysis

    [vec, val] = eig(M);
    val = diag(val);

    real = (imag(val) == 0);

    [~, i1] = max(abs(val));
    [~, i2] = min(abs(val));

    if(~real(i1) || ~real(i2))

        error('Imaginary eigenvalues are dominant')

    end

    % Stable and unstable eigenvectors

    a = 1:6;
    unstable = vec(:, a(i1));
    stable = vec(:, a(i2));

    % Preallocate arrays
    XS_in = zeros(6, N); 

    XS_out = XS_in; 
    XU_in = XS_in; 
    XU_out = XS_in;

    % Apply perturbations to each of N points

    for ii = 1:N

        % Grab state transition
        Phi = reshape(Y(ii, 7:42), 6, 6);

        % Grab state
        y = Y(ii, 1:6)';

        % Map stable and unstable vectors forward
        S = Phi*stable;
        U = Phi*unstable;

        % Get unit vectors; normalise

        S = S/norm(S);  U = U/norm(U);

        % Create perturbation vector

        pert = xpert;                                                           % Perturb by 10-6

        % Perturb conditions
        
        XS_in(:, ii) = y  + S.*pert;
        XS_out(:,ii) = y - S .* pert; 

    end

end
