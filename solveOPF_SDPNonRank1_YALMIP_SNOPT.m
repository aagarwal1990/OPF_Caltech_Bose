% solveOPF.m
%
% Solves OPF through:
%    SNOPT - with YAMLIP parser

% clear all
% close all
% clc

case_num = 'case39';
file_name = strcat(case_num, 'nonRank1_vars.mat');
use_line_limits = 1;
load(file_name);

[PgMax, PgMin, QgMax, QgMin, Pd, Qd, Fmax, conditionObj, costGen2, ...
 costGen1, costGen0, WMax, WMin, Phi, Psi, JJ, Ff, Tt, n, m] ...
 = setUpOptimVar(case_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SNOPT - SQP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting point for SQP
V_k = V0';
lambda_k = lambda0;
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

exp_V = sdpvar(2*n, 1); 
Pg    = sdpvar(n, 1); 
Qg    = sdpvar(n, 1);
Pinj  = sdpvar(n, 1); 
Qinj  = sdpvar(n, 1); 
Vsq   = sdpvar(n, 1);
aux   = sdpvar(n, 1);
Pf    = sdpvar(m, 1);
Pt    = sdpvar(m, 1);

OBJECTIVE    = sum(aux);
CONSTRAINTS = [];
for kk = 1:n
    CONSTRAINTS = [CONSTRAINTS, Pinj(kk) == exp_V'* exp_Phi{kk} * exp_V];
    CONSTRAINTS = [CONSTRAINTS, Qinj(kk) == exp_V'* exp_Psi{kk} * exp_V];
    CONSTRAINTS = [CONSTRAINTS, Vsq(kk)  == (exp_V(kk))^2 + (exp_V(kk + n))^2];
    CONSTRAINTS = [CONSTRAINTS, aux(kk)  >= costGen2(kk) * Pg(kk)^2 + + costGen1(kk) * Pg(kk) + costGen0(kk)];
end
        
% Line limits
for bb = 1:m
    CONSTRAINTS = [CONSTRAINTS, Pf(bb) == exp_V' * exp_Ff{bb} * exp_V];
    CONSTRAINTS = [CONSTRAINTS, Pt(bb) == exp_V' * exp_Tt{bb} * exp_V];
end

CONSTRAINTS = [CONSTRAINTS, Pinj == Pg - Pd];
CONSTRAINTS = [CONSTRAINTS, Qinj == Qg - Qd];

CONSTRAINTS = [CONSTRAINTS, Pg - PgMax <= 0];
CONSTRAINTS = [CONSTRAINTS, PgMin - Pg <= 0];
CONSTRAINTS = [CONSTRAINTS, Qg - QgMax <= 0];
CONSTRAINTS = [CONSTRAINTS, QgMin - Qg <= 0];
CONSTRAINTS = [CONSTRAINTS, Vsq - WMax <= 0];
CONSTRAINTS = [CONSTRAINTS, WMin - Vsq <= 0];

if use_line_limits == 1
    CONSTRAINTS = [CONSTRAINTS, Pf - line_limits <= 0];
    CONSTRAINTS = [CONSTRAINTS, -line_limits - Pt <= 0];
end

assign(exp_V, exp_V_k);
OPTIONS = sdpsettings('solver','snopt','usex0', 1);

% exp_V_random = ones(2*n, 1) * 1;
% assign(exp_V, exp_V_random);
% OPTIONS = sdpsettings('solver','snopt', 'usex0', 1);

%%%%%%%%%%% Solve using SNOPT %%%%%%%%%%%%%%%%%%%%%%%
solvesdp(CONSTRAINTS,OBJECTIVE, OPTIONS);
objective_value = sum(double(aux)) * conditionObj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1) Check if current V_k minimizes objective value 
% 2) Satisfies constraint set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_cand = complex(double(exp_V(1:n)), double(exp_V(n+1:2*n)));

% Set up necessary variables to check constraints 
epsilon = 10^-4;
assert(min(real(double(Pg)) - PgMax <= epsilon)==1)
assert(min(PgMin - real(double(Pg)) <= epsilon)==1)
assert(min(real(double(Qg)) - QgMax <= epsilon)==1)
assert(min(QgMin - real(double(Qg)) <= epsilon)==1)
assert(min(double(Vsq) - WMax <= epsilon)==1)
assert(min(WMin - double(Vsq) <= epsilon)==1)
if use_line_limits == 1
    assert(min(double(Pf) - line_limits <= epsilon)==1)
    assert(min(-line_limits - double(Pt) <= epsilon)==1)
end
sprintf('SNOPT - SQP: constraints satisfied with epsilon = %d', epsilon)
    
