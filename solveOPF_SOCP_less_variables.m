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

run_matpower = 0;
display('Check whether loadcase is commented')
%%%%%%%%%%%% COMMENT OUT FOR LOOP
% clear all
% close all
% clc
% 
% mpc = loadcase('case9');
%%%%%%%%%%%%

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
costGen1    = zeros(n, 1);
costGen0    = zeros(n, 1);

costGen2(genBuses) ...
            = mpc.gencost(:, 5) * (mpc.baseMVA^2) / conditionObj;
costGen1(genBuses) ...
            = mpc.gencost(:, 6) * mpc.baseMVA / conditionObj;
costGen0(genBuses) ...
            = mpc.gencost(:, 7) / conditionObj;
        
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
    JJ{k}       = (e_mat(:,k) * (e_mat(:, k)'));
end
clear y_k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SOCP relaxation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('--------- SOCP - less variables calculation ----------')

%start time
tic;

cvx_begin
    variables Pg(n) Qg(n) Pinj(n) Qinj(n) Vsq(n) aux(n); 
    variables Pf(m) Pt(m);
    variable W(n, n) hermitian
    minimize sum(aux)
    subject to
        
        
        for kk = 1:n
            var_1 = real(Phi{kk}(kk, kk) * W(kk, kk));
            var_2 = real(Psi{kk}(kk, kk) * W(kk, kk));
            
            for i = 1:length(mpc.branch(:, 1)) 
                if mpc.branch(i, 1) == kk
                    to = mpc.branch(i, 2); 
                    var_1 = var_1 + real(Phi{kk}(kk, to) * W(to, kk));
                    var_1 = var_1 + real(Phi{kk}(to, kk) * W(kk, to));
%                     display(mpc.branch(i, 1));
%                     display(mpc.branch(:, 1));
%                     display(mpc.branch(:, 2));
%                     display(Phi{kk}(kk, to));
%                     display(Phi{kk}(to, kk));
%                     display(Phi{kk});
                    var_2 = var_2 + real(Psi{kk}(kk, to) * W(to, kk));
                    var_2 = var_2 + real(Psi{kk}(to, kk) * W(kk, to));
                end
                if mpc.branch(i, 2) == kk
                    from = mpc.branch(i, 1);
                    var_1 = var_1 + real(Phi{kk}(from, kk) * W(kk, from));
                    var_1 = var_1 + real(Phi{kk}(kk, from) * W(from, kk));
%                     display(mpc.branch(i, 2));
%                     display(mpc.branch(:, 1));
%                     display(mpc.branch(:, 2));
%                     display(Phi{kk}(kk, to));
%                     display(Phi{kk}(to, kk));
%                     display(Phi{kk});
                    var_2 = var_2 + real(Psi{kk}(from, kk) * W(kk, from));
                    var_2 = var_2 + real(Psi{kk}(kk, from) * W(from, kk));
                end
            end
                                            
            Pinj(kk) == var_1;
            Qinj(kk) == var_2;
            Vsq(kk)  == W(kk, kk);

            costGen2(kk) * Pg(kk)^2 ...
                + costGen1(kk) * Pg(kk) ...
                + costGen0(kk) <= aux(kk);
            
        end
    
        Pinj == Pg - Pd;
        Qinj == Qg - Qd;
    
        
        Pg <= PgMax;
        Pg >= PgMin;
        Qg <= QgMax;
        Qg >= QgMin;
        Vsq >= WMin;
        Vsq <= WMax;
        
%         % Line limits
%         for bb = 1:m
%             Pf(bb) == real(trace(Ff{bb} * W));
%             Pt(bb) == real(trace(Tt{bb} * W));
%         end
%         
%                 
%         Pf <= Fmax;
%         Pt <= Fmax;


        for i = 1:length(mpc.branch(:, 1))
            fr = mpc.branch(i:i, 1);
            to = mpc.branch(i:i, 2);
            temp_matrix = [[W(fr, fr), W(fr, to)]', [W(to, fr), W(to, to)]'];
            temp_matrix == hermitian_semidefinite( 2 );
        end
         
cvx_end

toc
elapsed_time = toc;

display(strcat('Total cost = ', num2str(sum(aux)*conditionObj)));
display(' '); display(' '); display(' ');


eig_ratio_lst = zeros(1,m); % record max eigs in each matrix
display('Eigenvalue_1, Eigenvalue_2,       Eigenvalue_Ratio');
for i = 1:length(mpc.branch(:, 1))
    fr = mpc.branch(i:i, 1);
    to = mpc.branch(i:i, 2);
    mat = [[W(fr, fr), W(fr, to)]', [W(to, fr), W(to, to)]'];
    eig_lst = eigs(mat);
    disp([num2str(eig_lst(1)) blanks(10) num2str(eig_lst(2)) blanks(10) num2str(eig_lst(2)/eig_lst(1))]);
    eig_ratio_lst(ii) = eig_lst(2)/eig_lst(1);
end


%%%%%%%%%%% code added by karan
objective_value = sum(aux)*conditionObj;
maxEigRatio = max(eig_ratio_lst)

%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Matpower's solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_matpower == 1
    display('--------- MATPOWER optimization ----------')
    opt = mpoption('OPF_FLOW_LIM', 1);
    results = runopf(mpc, opt);

    obj = results.f;
    display(strcat('Total cost = ', num2str(obj)));
end

% % % % 
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
