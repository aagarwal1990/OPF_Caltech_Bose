% solveOPF.m
%
% Solves OPF through 2 methods:
%   1. SDP relaxation
%   2. Matpower's internal solver.
%
% Author: Subhonmesh Bose.
%
% Requires Matpower, CVX and SeDuMi.

clear all
close all
clc

case_num = 'case14';
%%%%%%%%%%%%

display('\n');
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
%         lam7 : Pf - line_limits <= 0;
%         lam8 : -line_limits - Pt <= 0;
        
        W == hermitian_semidefinite( n );
cvx_end

% get elapsed time
toc
elapsed_time = toc;

% get objective value
objective_value = sum(aux) * conditionObj

% get max eig ratio
eig_lst = eig(W);
max_eig = max(eig_lst);
maxEigRatio = max(eig_lst(eig_lst ~= max_eig))/max_eig

% get voltage values
[vec, lamda] = eigs(W);
eig_1 = lamda(1);
R = chol(W);
V0 = R(1, :);

sprintf('Checking assertions on SDP V0')
epsilon = 10^-4;
assert(min(real(Pg) - PgMax <= epsilon)==1)
assert(min(PgMin - real(Pg) <= epsilon)==1)
assert(min(real(Qg) - QgMax <= epsilon)==1)
assert(min(QgMin - real(Qg) <= epsilon)==1)
assert(min(Vsq - WMax <= epsilon)==1)
assert(min(WMin - Vsq <= epsilon)==1)
% % assert(min(Pf - line_limits <= epsilon)==1)
% % assert(min(-line_limits - Pt <= epsilon)==1)
% sprintf('constraints are satsified with epsilon = %d', epsilon)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Run Matlab FMINCON
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_k = V0';
exp_V_k = cat(1,real(V_k), imag(V_k));

options = optimset('GradObj','on','GradConstr','on', 'Display', 'iter', 'Algorithm','sqp');
[x,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon('objfun',exp_V_k,[],[],[],[],[],[],'constraints', options);

V_fin = complex(x(1:n), x(n+1:2*n))
objective_value = sum(fval) * conditionObj

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Manual Check of Feasibility Constraints
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lambda0 = [lam1', lam2', lam3', lam4', lam5', lam6']%, lam7', lam8'];
% 
% [hess_lagrangian, objective_value, V_fin]  = runSQP( V0, lambda0', case_num );
% objective_value
% V_fin

% Checking if V_k satisfies original constraints
display('Checking QCQP constraints with V_fin')
Pinj = zeros(n,1);
Qinj = zeros(n,1);
Vsq = zeros(n,1);
Pf = zeros(m,1);
Pt = zeros(m,1);
for kk = 1:n
    Pinj(kk) = V_fin' * Phi{kk} * V_fin;
    Qinj(kk) = V_fin' * Psi{kk} * V_fin;
    Vsq(kk)  = abs(V_fin(kk));
end

Pg = Pinj + Pd;
Qg = Qinj + Qd;

% Line limits
for bb = 1:m
    Pf(bb) = V_fin' * Ff{bb}* V_fin;
    Pt(bb) = V_fin' * Tt{bb}* V_fin;
end

% Contraints
epsilon = 10^-4;
assert(min(real(Pg) - PgMax <= epsilon)==1)
assert(min(PgMin - real(Pg) <= epsilon)==1)
assert(min(real(Qg) - QgMax <= epsilon)==1)
assert(min(QgMin - real(Qg) <= epsilon)==1)
assert(min(Vsq - WMax <= epsilon)==1)
assert(min(WMin - Vsq <= epsilon)==1)
% assert(min(Pf - line_limits <= epsilon)==1)
% assert(min(-line_limits - Pt <= epsilon)==1)
sprintf('constraints satisfied with epsilon = %d', epsilon)