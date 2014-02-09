function [hess_lagrangian, objective_value, V_fin ] = runSQP( V0, lambda0, case_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% solveOPF.m
%
% Solves OPF through 2 methods:
%   1. SDP relaxation
%   2. Matpower's internal solver.
%
% Author: Subhonmesh Bose.
%
% Requires Matpower, CVX and SeDuMi.


mpc = loadcase(case_num);
n           = size(mpc.bus, 1);
m           = size(mpc.branch, 1);
genBuses    = mpc.gen(:, 1);
ng          = size(genBuses, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Change bus order
%%%%%%%%%%%%%%%%%%%%%%%%%


mpc.bus         = sortrows(mpc.bus);
mapBus          = mpc.bus(:, 1);
mpc.bus(:, 1)   = (1:n)';
for jj = 1:m
    mpc.branch(jj, 1) = find(mapBus == mpc.branch(jj, 1));
    mpc.branch(jj, 2) = find(mapBus == mpc.branch(jj, 2));
end
for jj = 1:ng
    mpc.gen(jj, 1) = find(mapBus == mpc.gen(jj, 1));
end

clear mapBus jj 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the data up for generation and demand. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PgMax       = zeros(n, 1);
PgMin       = zeros(n, 1);

QgMax       = zeros(n, 1);
QgMin       = zeros(n, 1);

PgMax(genBuses) ...
            = mpc.gen(:, 9) / mpc.baseMVA;
PgMin(genBuses) ...
            = mpc.gen(:, 10) / mpc.baseMVA;
QgMax(genBuses) ...
            = mpc.gen(:, 4) / mpc.baseMVA;
QgMin(genBuses) ...
            = mpc.gen(:, 5) / mpc.baseMVA;
    
Pd          = mpc.bus(:, 3) / mpc.baseMVA;
Qd          = mpc.bus(:, 4) / mpc.baseMVA;

Fmax        = mpc.branch(:, 6) / mpc.baseMVA;


conditionObj = 10 * mpc.baseMVA;

costGen2    = zeros(n, 1);
costGen1    = ones(n, 1) * 3;
costGen0    = zeros(n, 1);

% costGen2(genBuses) ...
%             = mpc.gencost(:, 5) * (mpc.baseMVA^2) / conditionObj;
% costGen1(genBuses) ...
%             = mpc.gencost(:, 6) * mpc.baseMVA / conditionObj;
% costGen0(genBuses) ...
%             = mpc.gencost(:, 7) / conditionObj;
        
WMax        = mpc.bus(:, 12) .^ 2;
WMin        = mpc.bus(:, 13) .^ 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optimization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = sqrt(-1);

for bb = 1 : m
   if mpc.branch(bb, 3) == 0
       mpc.branch(bb, 3) = 10^(-4);           
   end                                    
end
clear bb

[Ybus, Yf, Yt] ...
            = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);

        
e_mat = eye(n);

Phi = {};
Psi = {};
JJ  = {};

% Cmatrices = {};
% bScalars  = {};

for k=1:n
    y_k         = e_mat(:,k) * (e_mat(:, k)') * Ybus;
    Phi{k}      = (ctranspose(y_k) + y_k)/2;
    Psi{k}      = (ctranspose(y_k) - y_k)/(2*j);
    JJ{k}       = e_mat(:,k) * (e_mat(:, k)');
end
clear y_k

Ff      = {};
Tt      = {};

for bb = 1:m

    eb      = zeros(m, 1);
    eb(bb)  = 1;
    
    eff     = zeros(n, 1);
    eff(mpc.branch(bb, 1)) ...
            = 1;
    
    ett     = zeros(n, 1);
    ett(mpc.branch(bb, 2)) ...
            = 1;
      
    Ff{bb}  = ctranspose(Yf) * eb * (eff');
    Ff{bb}  = (Ff{bb} + ctranspose(Ff{bb})) / 2;
    Tt{bb}  = ctranspose(Yt) * eb * (ett');
    Tt{bb}  = (Tt{bb} + ctranspose(Tt{bb})) / 2;

end
clear bb eb eff ett

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SQP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

line_limits = ones(m, 1) * 100;
% Impose line limit on Branch (7,8)
line_limits(14) = -0.995;
% Impose line limit on Branch (7,9)
line_limits(15) = 0.5000;

% Starting point for SQP
V_k = V0';
iter_diff = 100;
lambda_k = lambda0;
count = 0;
exp_V_k = cat(1,real(V_k), imag(V_k));
exp_Phi = {};
exp_Psi = {};
exp_Ff  = {};
exp_Tt  = {};
for kk = 1:n
    top         = cat(2, real(Phi{kk}), -imag(Phi{kk}));
    bottom      = cat(2, imag(Phi{kk}),  real(Phi{kk}));
    exp_Phi{kk} = cat(1, top, bottom);

    top         = cat(2, real(Psi{kk}), -imag(Psi{kk}));
    bottom      = cat(2, imag(Psi{kk}),  real(Psi{kk}));
    exp_Psi{kk} = cat(1, top, bottom);
end

for kk = 1:m
    top        = cat(2, real(Ff{kk}), -imag(Ff{kk}));
    bottom     = cat(2, imag(Ff{kk}),  real(Ff{kk}));
    exp_Ff{kk} = cat(1, top, bottom);
    
    top        = cat(2, real(Tt{kk}), -imag(Tt{kk}));
    bottom     = cat(2, imag(Tt{kk}),  real(Tt{kk}));
    exp_Tt{kk} = cat(1, top, bottom);
end

while and(iter_diff > 10^-4, count < 10)
    count = count + 1
    grad_cost = zeros(2*n,1);
    for kk = 1:n
        grad_cost = grad_cost + 2*costGen1(kk)*(exp_Phi{kk}*exp_V_k);
    end
    
    % Calculation of jacobian of lagrangian
    jacobian_g = zeros(6*n + 2*m, 2*n);
    
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
    for kk = 1:m
        jacobian_g(6*n+kk,:)       =  2*(exp_Ff{kk}*exp_V_k)';
        jacobian_g(6*n + m + kk,:) = -2*(exp_Tt{kk}*exp_V_k)';
    end
    
    grad_lagrangian = grad_cost' + lambda_k' * jacobian_g;
    
    % Calculation of hessian of lagrangian
    hess_lagrangian = zeros( 2*n,2*n);
    
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
    for kk = 1:m
         hess_lagrangian = hess_lagrangian ...
                          + 2*lambda_k(6*n+kk, :)   * exp_Ff{kk} ...
                          - 2*lambda_k(6*n+m+kk, :) * exp_Tt{kk};
    end
    
    cvx_begin quiet
        variables exp_V(2*n) V(n) W(n, n) obj;
        variables Pg(n) Qg(n) Pinj(n) Qinj(n) Vsq(n) aux(n); 
        variables Pf(m) Pt(m);
        dual variables lam1 lam2 lam3 lam4 lam5 lam6 lam7 lam8;
        minimise obj;
        subject to
            obj >= grad_lagrangian * (exp_V-exp_V_k) ...
                        + 1/2 * (exp_V-exp_V_k)' * hess_lagrangian ...
                        * (exp_V - exp_V_k);
            for kk = 1:n
                Pinj(kk) == exp_V_k'* exp_Phi{kk} * exp_V_k ;
                Qinj(kk) == exp_V_k'* exp_Psi{kk} * exp_V_k;
                Vsq(kk)  == (exp_V_k(kk))^2 + (exp_V_k(kk + n))^2;
            end

            Pinj == Pg - Pd;
            Qinj == Qg - Qd;
            
            % Line limits
            for bb = 1:m
                Pf(bb) == exp_V_k' * exp_Ff{bb} * exp_V_k;
                Pt(bb) == exp_V_k' * exp_Tt{bb} * exp_V_k;
            end
            
            % Contraints
            lam1 : Pg - PgMax + jacobian_g(1:n, :)      *(exp_V-exp_V_k)<= 0;
            lam2 : PgMin - Pg + jacobian_g(n+1:2*n, :)  *(exp_V-exp_V_k)<= 0;
            lam3 : Qg - QgMax + jacobian_g(2*n+1:3*n, :)*(exp_V-exp_V_k)<= 0;
            lam4 : QgMin - Qg + jacobian_g(3*n+1:4*n, :)*(exp_V-exp_V_k)<= 0;
            lam5 : Vsq - WMax + jacobian_g(4*n+1:5*n, :)*(exp_V-exp_V_k)<= 0;
            lam6 : WMin - Vsq + jacobian_g(5*n+1:6*n, :)*(exp_V-exp_V_k)<= 0;
            lam7 : Pf - line_limits + ...
                                jacobian_g(6*n+1:6*n+m, :)*(exp_V-exp_V_k)<= 0;
            lam8 : - line_limits - Pt + ...
                                jacobian_g(6*n+1+m:6*n+2*m, :)*(exp_V-exp_V_k)<= 0;

    cvx_end
    
    iter_diff = norm(exp_V - exp_V_k);
    
    exp_V_k = exp_V;
    
    lambda_temp = [lam1', lam2', lam3', lam4', lam5', lam6', lam7', lam8'];
    lambda_k = lambda_temp';
    
end

objective_value = zeros(n, 1);
for kk = 1:n
   objective_value(kk) = costGen2(kk) * Pg(kk)^2 ...
        + costGen1(kk) * Pg(kk) ...
        + costGen0(kk);
end
V_fin = complex(exp_V_k(1:n), exp_V_k(n+1:2*n));
objective_value = sum(objective_value)*conditionObj;
end