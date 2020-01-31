% Equations of Motion for the three-dimensional CR3BP
%
% List of inputs
% --------------
% x0: Initial state vector
% mu: System mass parameter
% 
% Outputs
% -------
% x, xdot: State vector and state vector derivative
%x

function [out1, out2, out3] = threeBodyEOM_discretize(~, x, flag)

    if nargin < 3 || isempty(flag)

        mu = 3.0032080443e-06; 
        x0 = x(1:6, 1);
        xdot = zeros(7,1);
        r1 = sqrt((x0(1)+mu)^2+x0(2)^2+x0(3)^2);
        r2 = sqrt((x0(1)-1+mu)^2+x0(2)^2+x0(3)^2);
        xdot(1:3) = x0(4:6);

        % Now, variational equations

        C1 = (1-mu)/r1^3;
        C2 = mu/r2^3;
        xdot(4) = 2*x0(5)+x0(1) - C1*(x0(1)+mu) - C2*(x0(1) - 1 + mu);
        xdot(5) = -2*x0(4) + x0(2) - C1*x0(2) - C2*x0(2);
        xdot(6) = -C1*x0(3) - C2*x0(3);
        out1 = xdot;
        

    else

        switch flag

            case 'events'
                
                    % Value - array of zeros for now
                    
                    if (x(7) == -1)
                        
                        value = x(2);
                    
                    else
                        
                        value = x(7) - abs(atan2(x(2), x(1)));
                        
                    end
                       
                    isterminal = 1;
                    direction = 0;
                    
                    out1 = value;
                    out2 = isterminal;
                    out3 = direction;

        end

    end

end