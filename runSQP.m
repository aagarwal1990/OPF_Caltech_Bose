function [ objective_value, V_fin ] = runSQP( V0, lamda0, case_num )
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
V_k = V0;
iter_diff = 100;
lamda_k = lamda0;
count = 0
exp_V_k = cat(1,real(V_k), imag(V_k));
exp_Phi = {};
exp_Psi = {};
exp_Ff  = {};
exp_Tt  = {};
for kk = 1:n
    top = cat(2, real(Phi{kk}), -imag(Phi{kk}));
    bottom = cat(2, imag(Phi{kk}), real(Phi{kk}));
    exp_Phi{kk} = cat(1, top, bottom);

    top = cat(2, real(Psi{kk}), -imag(Psi{kk}));
    bottom = cat(2, imag(Psi{kk}), real(Psi{kk}));
    exp_Psi{kk} = cat(1, top, bottom);
end

for kk = 1:m
    top = cat(2, real(Ff{kk}), -imag(Ff{kk}));
    bottom = cat(2, imag(Ff{kk}), real(Ff{kk}));
    exp_Ff{kk} = cat(1, top, bottom);
    
    top = cat(2, real(Tt{kk}), -imag(Tt{kk}));
    bottom = cat(2, imag(Tt{kk}), real(Tt{kk}));
    exp_Tt{kk} = cat(1, top, bottom);
end

while and(iter_diff > 10^-4, count < 10)
    count = count + 1
    grad_cost = zeros(2*n,1);
    for kk = 1:n
        grad_cost = grad_cost + 2*costGen1(kk)*(exp_Phi{kk}*exp_V_k);
    end

    jacobian_g = zeros(6*n + 2*m, 2*n);

    % adding determinant of Pg <= PgMax

    for kk = 1:n
        jacobian_g(kk,:)     = 2*(exp_Phi{kk}*exp_V_k)';
        jacobian_g(n+kk,:)   = -2*(exp_Phi{kk}*exp_V_k)';
        jacobian_g(2*n+kk,:) = 2*(exp_Psi{kk}*exp_V_k)';
        jacobian_g(3*n+kk,:) = -2*(exp_Psi{kk}*exp_V_k)';
        jacobian_g(4*n+kk,:) = 2*exp_V_k';
        jacobian_g(5*n+kk,:) = -2*exp_V_k';
    end
    
    for kk = 1:m
        jacobian_g(6*n+kk,:) = 2*(exp_Ff{kk}*exp_V_k)';
        jacobian_g(6*n + m + kk,:) = -2*(exp_Tt{kk}*exp_V_k)';
    end
    
    grad_lagrangian = grad_cost + lambda_k'*jacobian_g;

    hess_lagrangian = zeros(2*n,2*n);
    
    for kk = 1:n
        hess_lagrangian = hess_lagrangian + 2*costGen1(kk)*exp_Phi{kk} 
    end
    
    for kk = 1:m
        hess_lagrangian = hess_lagrangian + 2*(exp_Ff{kk} - exp_Tt{kk});
    end

    cvx_begin
        variable exp_V(2*n) V(n) W(n, n) obj;
        variables Pg(n) Qg(n) Pinj(n) Qinj(n) Vsq(n) aux(n); 
        variables Pf(m) Pt(m);
        dual variables lam1 lam2 lam3 lam4 lam5 lam6 lam7 lam8;
        minimise obj;
        subject to

            obj = grad_lagrangian'* (exp_V-exp_V_k) + 1/2 * (exp_V-exp_V_k)'* hess_lagrangian * (exp_V-exp_V_k);
            
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
            lam1 : Pg - PgMax <= 0;
            lam2 : PgMin - Pg <= 0;
            lam3 : Qg - QgMax <= 0;
            lam4 : QgMin - Qg <= 0;
            lam5 : Vsq - WMax <= 0;
            lam6 : WMin - Vsq <= 0;
            lam7 : Pf - line_limits <= 0;
            lam8 : -line_limits - Pt <= 0;

    cvx_end
    
    iter_diff = norm(exp_V - exp_V_k);
    exp_V_k = exp_V;
    lamda_temp = [lam1', lam2', lam3', lam4', lam5', lam6', lam7', lam8'];
    for i = 1:6
        index = n * i;
        lamda_k(index - n + 1 : index, 1) = lamda_temp(i);
    end

    for i = 7:8
        index = m * i;
        lamda_k(index - m + 1 : index, 1) = lamda_temp(i);
    end
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




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new code to determine sj
% from_verts = mpc.branch(:,1);
% to_verts = mpc.branch(:,2);
% 
% res = mpc.branch(:,3);
% rctn = mpc.branch(:,4);
% 
% 
% Adj = sparse(n,n);
% Admittance = sparse(n,n)
% for ii = 1:m
%     jj = from_verts(ii);
%     kk = to_verts(ii);
%     Adj(jj, kk) = 1;
%     Adj(kk, jj) = 1;
%     
%     y_jk = complex(res(ii), -rctn(ii))/sqrt(res(ii)^2 + rctn(ii)^2)
%     Admittance(jj, kk) = y_jk;
%     Admittance(kk, jj) = y_jk;
% end
% 
% % to ensure that diagonal entries are zero
% Adj(logical(eye(n)))=0;
% 
% rAdmittance = real(Admittance);% real part of admittance matrix
% iAdmittance = imag(Admittance);% imaginary part of admittance matrix
% 
% nbrs = cell(1,n);
% for jj = 1:n
%     nbrs{jj} = find(Adj(jj,:));
% end
% 
% 
% % Following is code for iterative step which takes as input V_k
% % iterative step to determine s for each iteration
% % assuming V_k (column vector) has been determined
% 
% 
% 
% % % matrix containing jth row multiplied by (V_k)_j
% % A = Adj * diag(V_k);
% % B = ((temp-temp').*Admittance) * diag(V_k);
% % s = B(logical(eye(n))); % this a part of g(x_k) in the book chapter 4
% % 
% % num_ineq = 2*(n+m); % min max for voltages and branch powers
% % 
% % lambda = zeros(num_ineq, 1);
% 
% 
% 
% 
% % To determine the cost function which is dependent only on the real part
% % as is seen in the solveOPF cvx, since cost is function of Pg only
% 
% s = zeros(n,1);
% 
% for jj = 1:n
%     s(jj) = abs(V_k(jj))^2*sum(Admittance(:,jj)) - V_k(jj)*sum(Admittance(jj,:)*V_k);
% end
% Pinj = real(s);
% Qinj = imag(s);
% 
% Vsq = abs(V_k).^2;
% rV_k = real(V_k);% real part of voltage vector
% iV_k = imag(V_k);% imaginary part of voltage vector
% % no we have the constant parts of the inequalities . Need to evaluate the linear parts
% % for that we use jacobian matrix in the various components
% 
% % To determine the jacobian matrices
% jacobian_Pinj = sparse(2*n, n);
% jacobian_Qinj = sparse(2*n, n);
% for jj = 1:n
%     
%     jacobian_Pinj(2*jj-1,jj) = 2*rV_k(jj)*sum(rAdmittance(jj,:)) ...
%         + rAdmittance(jj,:)*rV_k ...
%         - iV_k(jj)*sum(iAdmittance(jj,:)) ...
%         + iAdmittance(jj,:)*iV_k ...
%         + iV_k(jj)*rAdmittance(JJ,:)*iV_k; % d /d real()
%     jacobian_Pinj(2*jj, jj) = -2*iV_k(jj)*rAdmittance(jj,:)*iV_k ...
%         + rV_k(jj)*iAdmittance(jj,:)*rV_k ...
%         - rV_k(jj)*iAdmittance(jj:1)*i
%     jacobian_Qinj
%     jacobian_Qinj
%     
%     nbrs_jj = nbrs{jj}
%     for ii = 1:size(nbrs_jj, 2)
%         jacobian_Pinj(2*jj-1, ii) = 
%         jacobian_Pinj(2*jj, ii) = 
%         jacobian_Qinj
%         jacobian_Qinj
%     end
% end