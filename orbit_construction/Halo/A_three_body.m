function [A, stable, unstable] = A_three_body(X, mu)
    
    % Computes the Jacobian of the state derivative using the 
    % state vector X, and return the stable and unstable eigenvectors.

    % Ensure X was passed in as a column vector
    if(size(X,2) > 1), X = X'; end

    % Vectors between primaries and satellite
    r1 = X(1:3);
    r1(1) = r1(1) + mu;

    r2 = r1;
    r2(1) = r2(1) - 1;

    % Vector norms
    R1 = norm(r1);
    R2 = norm(r2);

    % Form second derivatives matrices
    B = diag([1 1 0]) - ((1-mu)/R1^3 + mu/R2^3)*eye(3) + ...
        3*(1-mu)*(r1*r1')/R1^5 + 3*mu*(r2*r2')/R2^5;

    C = [0 2 0; -2 0 0; 0 0 0];

    % Form entire A matrix
    A = [zeros(3), eye(3); B, C];

    % Eigenvalues/vectors
    if(nargout > 1)
        % Gather the six eigenvalues and vectors
        [vecs, vals] = eig(A);
        
        % Look for the positive and negative real eigenvalues
        vals = diag(vals);
        pos_ind = (imag(vals) == 0) & (vals > 0);
        neg_ind = (imag(vals) == 0) & (vals < 0);
        
        if(sum(pos_ind) == 1 && sum(neg_ind) == 1)
            % Convert to indices and extract eigenvectors
            a = 1:6;
            pos = a(pos_ind);
            neg = a(neg_ind);
            
            unstable = vecs(:, pos);
            stable = vecs(:, neg);
            
        else
            disp(vals)
            error('More than two real eigenvalues - check form of A matrix') 
        end
    end

    end