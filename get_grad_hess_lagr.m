function [jacobian_g, grad_lagrangian, hess_lagrangian] = ...
    get_grad_hess_lagr(costGen0, costGen1, costGen2, ...
                       exp_Phi, exp_Psi, exp_Ff, exp_Tt, exp_V_k, ...
                       lambda_k, n, m, use_line_limits)
                   
grad_cost = zeros(2*n,1);
for kk = 1:n
    grad_cost = grad_cost + 2*costGen1(kk)*(exp_Phi{kk}*exp_V_k);
end

% Calculation of jacobian of lagrangian
if use_line_limits == 1
    jacobian_g = zeros(6*n + 2*m, 2*n);
else
    jacobian_g = zeros(6*n, 2*n);
end

% bus constraints segment of jacobian of lagrangian
for kk = 1:n
    jacobian_g(kk,:)     =  2*(exp_Phi{kk}*exp_V_k)';
    jacobian_g(n+kk,:)   = -2*(exp_Phi{kk}*exp_V_k)';
    jacobian_g(2*n+kk,:) =  2*(exp_Psi{kk}*exp_V_k)';
    jacobian_g(3*n+kk,:) = -2*(exp_Psi{kk}*exp_V_k)';
    jacobian_g(4*n+kk,:) =  2*exp_V_k';
    jacobian_g(5*n+kk,:) = -2*exp_V_k';
end

% branch constraints segment of jacobian of lagrangian
if use_line_limits == 1
    for kk = 1:m
        jacobian_g(6*n+kk,:)       =  2*(exp_Ff{kk}*exp_V_k)';
        jacobian_g(6*n + m + kk,:) = -2*(exp_Tt{kk}*exp_V_k)';
    end
end

grad_lagrangian = grad_cost' + lambda_k' * jacobian_g;

% Calculation of hessian of lagrangian
hess_lagrangian = zeros(2*n, 2*n);

% generation cost segment of lagrangian
for kk = 1:n
    hess_lagrangian = hess_lagrangian ...
            + 2 * costGen1(kk) * exp_Phi{kk};
end

% bus constraints segment of hessian of lagrangian
for kk = 1:n
    hess_lagrangian = hess_lagrangian ...
                        + 2 * lambda_k(kk, :)     * exp_Phi{kk} ...
                        - 2 * lambda_k(n+kk, :)   * exp_Phi{kk} ...
                        + 2 * lambda_k(2*n+kk, :) * exp_Psi{kk} ...
                        - 2 * lambda_k(3*n+kk, :) * exp_Psi{kk} ...
                        + 2 * lambda_k(4*n+kk, :) * eye(2*n)    ...
                        - 2 * lambda_k(5*n+kk, :) * eye(2*n);
end

% branch constraints segment of hessian of lagrangian
if use_line_limits == 1
    for kk = 1:m
         hess_lagrangian = hess_lagrangian ...
                          + 2*lambda_k(6*n+kk, :)   * exp_Ff{kk} ...
                          - 2*lambda_k(6*n+m+kk, :) * exp_Tt{kk};
    end
end