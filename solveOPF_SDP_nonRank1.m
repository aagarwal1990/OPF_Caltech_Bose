% solveOPF.m
%
% Solves OPF through 2 methods:
%   1. SDP relaxation
%   2. Matpower's internal solver.
%
% Author: Subhonmesh Bose.
%
% Requires Matpower, CVX and SeDuMi.
%

display('Check whether loadcase is commented')
%%%%%%%%%%%%   COMMENT OUT FOR LOOP
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

cvx_begin
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
        lam7 : Pf - line_limits <= 0;
        lam8 : -line_limits - Pt <= 0;
        
        W == hermitian_semidefinite( n );
cvx_end

% get elapsed time
toc
elapsed_time = toc;

% get objective value
objective_value = sum(aux)*conditionObj;

% get max eig ratio
eig_lst = eig(W);
max_eig = max(eig_lst);
maxEigRatio = max(eig_lst(eig_lst ~= max_eig))/max_eig;
eigs(W)

% get voltage values
[vec, lamda] = eigs(W);
eig_1 = lamda(1);
R = chol(W);
V0 = R(1, :);
lamda0 = {lam1, lam2, lam3, lam4, lam5, lam6, lam7, lam8};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Matpower's solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if run_matpower == 1
%     display('--------- MATPOWER optimization ----------')
%     opt = mpoption('OPF_FLOW_LIM', 1);
%     results = runopf(mpc, opt);
% 
%     obj = results.f;
%     display(strcat('Total cost = ', num2str(obj)));
% end
% uA = 7.5/100;
% uB = 10/100;
% uC = 20/100;
% 20000*(5.05 + uA+uB+2*uC)
% 
% u = [uA, uB, uC]';
% k = ones(3);
% 
% varA = (7/100)^2;
% varB = (12/100)^2;
% varC = (18/100)^2;
% 
% corelAB = 0.7;
% corelCA = -0.5;
% corelBC = -0.3;
% 
% covarAB = corelAB*sqrt(varA)*sqrt(varB);
% covarBC = corelBC*sqrt(varC)*sqrt(varB);
% covarCA = corelCA*sqrt(varC)*sqrt(varA);
% 
% tot_var = varA+varB+4*varC+2*covarAB+4*(covarCA+covarBC);
% 20000*tot_var
% sqrt(20000*tot_var)
% 
% covar_mat = [varA, covarAB, covarCA; covarAB, varB, covarBC; covarCA, covarBC, varC]
% covar_mat2 = cat(2,covar_mat,[0,0,0]');
% covar_mat2 = cat(1,covar_mat2,[0,0,0,0]);

% % % % if strcmp(cvx_status, 'Solved') ~= 1
% % % %     display('Problems in optimization');
% % % % end
% % % % 
% % % % % display('P power')
% % % % % display([Pg, PgMax])
% % % % % display('Q power')
% % % % % display([Qg, 0.6*Pg])
% % % % % display('Line flow')
% % % % % display([-Fmax, Pt, Pf, Fmax])
% % % % % display('Voltage')
% % % % % display([WMin, Vsq, WMax])
% % % % 
% % % % 
% % % % 
% % % % %display([mpc.branch(:,1), mpc.branch(:, 2), Pf, Pt, Fmax])
% % % % %Pd
% % % % %Pg
% % % % NF = sum(Pd);
% % % % 
% % % % if ~isnan(W)
% % % %     lambdas = sort(real(eig(W)), 1, 'descend');
% % % %     lambda123 = lambdas(1:3);
% % % % else
% % % %     lambda123 = NaN;
% % % % end
