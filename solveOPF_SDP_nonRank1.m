% solveOPF.m
%
% Solves OPF through 2 methods:
%   1. SDP relaxation
%   2. Matpower's internal solver.
%
% Author: Subhonmesh Bose.
%
% Requires Matpower, CVX and SeDuMi.

% clear all
% close all
% clc

case_num = 'case14';
use_line_limits =  1;
[PgMax, PgMin, QgMax, QgMin, Pd, Qd, Fmax, conditionObj, costGen2, ...
 costGen1, costGen0, WMax, WMin, Phi, Psi, JJ, Ff, Tt, n, m, bus, branch] ...
= setUpOptimVar(case_num);
        
line_limits = ones(m, 1) * 100;
% Impose line limit on Branch (7,8)
line_limits(14) = -0.995;
% Impose line limit on Branch (7,9)
line_limits(15) = 0.5000;
% start time
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SDP relaxation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display('--------- SDP calculation ----------')

cvx_begin quiet
    variables Pg(n) Qg(n) Pinj(n) Qinj(n) Vsq(n) aux(n); 
    variables Pf(m) Pt(m);
    dual variables lam1 lam2 lam3 lam4 lam5 lam6 lam7 lam8;
    variable W(n, n) hermitian
    minimize sum(aux)
    cvx_solver sedumi
    subject to
        
        
        for kk = 1:n
            Pinj(kk) == real( trace( Phi{kk} * W ));
            Qinj(kk) == real( trace( Psi{kk} * W ));
            Vsq(kk)  == W(kk, kk);
            
            costGen2(kk) * Pg(kk)^2 ...
                + costGen1(kk) * Pg(kk) ...
                + costGen0(kk) <= aux(kk);
        end
    
        Pinj == Pg - Pd;
        Qinj == Qg - Qd;
                
        % Line limits
        for bb = 1:m
            Pf(bb) == real(trace(Ff{bb} * W));
            Pt(bb) == real(trace(Tt{bb} * W));
        end
        
        % Contraints
        lam1 : Pg - PgMax <= 0;
        lam2 : PgMin - Pg <= 0;
        lam3 : Qg - QgMax <= 0;
        lam4 : QgMin - Qg <= 0;
        lam5 : Vsq - WMax <= 0;
        lam6 : WMin - Vsq <= 0;
        if use_line_limits == 1
            lam7 : Pf -line_limits  <= 0;
            lam8 : -Pt -line_limits <= 0;
        end
        
        W == hermitian_semidefinite( n );
cvx_end

% get elapsed time
toc
elapsed_time = toc;

% get objective value
objective_value_SDP = sum(aux) * conditionObj

% get max eig ratio
eig_lst = eig(W);
max_eig = max(eig_lst);
maxEigRatio = max(eig_lst(eig_lst ~= max_eig))/max_eig
line_summary = zeros(m, 3);
line_summary(:, 1) = Pf;
line_summary(:, 2) = branch(:, 1);
line_summary(:, 3) = branch(:, 2);
line_summary(:, 4) = 1:m;
line_summary;

% get voltage values
[vec, lamda] = eigs(W);
eig_1 = lamda(1);
R = chol(W);
V0 = R(1, :);

if use_line_limits == 1
    lambda0 = [lam1', lam2', lam3', lam4', lam5', lam6', lam7', lam8'];
else
    lambda0 = [lam1', lam2', lam3', lam4', lam5', lam6'];
end

file_name = strcat(case_num, 'nonRank1_vars.mat');
save(file_name, 'lambda0', 'V0', 'W', 'line_limits', 'objective_value_SDP');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Run SQP relaxation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% display('%%%%%% Starting SQP %%%%%%%%%%')
% [hess_lagrangian, objective_value, V_fin, num_iter]  = ...
%          runSQP( V0, lambda0', case_num, use_line_limits, line_limits);
% toc
% 
% objective_value
% num_iter;
% V_fin;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check if Constraints Satisfied
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%% SDP %%%%%%%%%%%%%%%%%%
% epsilon = 10^-4;
% assert(min(real(Pg) - PgMax <= epsilon)==1)
% assert(min(PgMin - real(Pg) <= epsilon)==1)
% assert(min(real(Qg) - QgMax <= epsilon)==1)
% assert(min(QgMin - real(Qg) <= epsilon)==1)
% assert(min(Vsq - WMax <= epsilon)==1)
% assert(min(WMin - Vsq <= epsilon)==1)
% if use_line_limits == 1
%     assert(min(Pf - line_limits <= epsilon)==1)
%     assert(min(-line_limits - Pt <= epsilon)==1)
% end
% sprintf('SDP: constraints are satsified with epsilon = %d', epsilon)
% 
% %%%%%%%%%%% SQP %%%%%%%%%%%%%%%%%%
% Pinj = zeros(n,1);
% Qinj = zeros(n,1);
% Vsq = zeros(n,1);
% Pf = zeros(m,1);
% Pt = zeros(m,1);
% for kk = 1:n
%     Pinj(kk) = V_fin' * Phi{kk} * V_fin;
%     Qinj(kk) = V_fin' * Psi{kk} * V_fin;
%     Vsq(kk)  = abs(V_fin(kk));
% 
% end
% 
% Pg = Pinj + Pd;
% Qg = Qinj + Qd;
% 
% % Line limits
% for bb = 1:m
%     Pf(bb) = V_fin' * Ff{bb}* V_fin;
%     Pt(bb) = V_fin' * Tt{bb}* V_fin;
% end
% 
% % Contraints
% epsilon = 10^-4;
% assert(min(real(Pg) - PgMax <= epsilon)==1)
% assert(min(PgMin - real(Pg) <= epsilon)==1)
% assert(min(real(Qg) - QgMax <= epsilon)==1)
% assert(min(QgMin - real(Qg) <= epsilon)==1)
% assert(min(Vsq - WMax <= epsilon)==1)
% assert(min(WMin - Vsq <= epsilon)==1)
% if use_line_limits == 1
%     assert(min(Pf - line_limits <= epsilon)==1)
%     assert(min(-line_limits - Pt <= epsilon)==1)
% end
% sprintf('SQP: constraints satisfied with epsilon = %d', epsilon)
