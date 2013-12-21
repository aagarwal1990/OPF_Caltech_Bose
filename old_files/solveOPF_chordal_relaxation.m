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


clear all
close all
clc

case_num = 'case30';
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

% Start Timer
tic

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
    Phi{k}      = sparse(ctranspose(y_k) + y_k)/2;
    Psi{k}      = sparse(ctranspose(y_k) - y_k)/(2*j);
    JJ{k}       = sparse(e_mat(:,k) * (e_mat(:, k)'));
end
clear y_k


% Get Maximal cliques for chordal extension
% Get neighbours 
[neighbours, M] = get_chordal_matrix_chole(mpc.branch(:, 1), mpc.branch(:, 2), n);

% Create Perfect Elimination Ordering of chordal extension
numbered_lst = false(1, n);
label_lst = zeros(1, n);
max_unNum_vertex = [0; 1];
PEO = zeros(1, n);

for i = n:-1:1
    node = max_unNum_vertex(2, 1);
    PEO(i) = node;  
    numbered_lst(node) = true;
    
    for j = 1:length(neighbours{node})
        neigh = neighbours{node}(j);
        label_lst(neigh) = label_lst(neigh) + 1;
    end
    
    max_label = 0;
    unNum_lst = find(~numbered_lst);
    for k = 1: length(unNum_lst);
        cand_node = unNum_lst(k);
        if label_lst(cand_node) > max_label
            max_label = label_lst(cand_node);
            max_unNum_vertex = [label_lst(cand_node); cand_node];
        end
    end
end

% Get Max Cliques
size = zeros(1, n);
max_clique = cell(1, n);
clique_index = 1;
M = {};

for i = 1:n
    node = PEO(i);
    if (isempty(neighbours{node}) == 1)
        max_clique{clique_index} = [node];
        clique_index = clique_index + 1;
    end
    
    X = {};
    for j = 1:length(neighbours{node})
        neigh = neighbours{node}(j);
        min_neigh_index = n;
        neigh_index = find(PEO == neigh);
        node_index = find(PEO == node);
        
        if neigh_index > node_index
            X{length(X) + 1} = neigh;
            if neigh_index < min_neigh_index
                min_neigh = neigh;
                min_neigh_index = neigh_index;
            end
        end
    end
    if (isempty(X) == 0)
        if size(node) < length(X)
            max_clique{clique_index} = [X{:}, node];
            clique_index = clique_index + 1;
        end 
        size(min_neigh) = max(size(min_neigh), length(X) - 1);
    end
end

max_clique = max_clique(1: clique_index - 1);
for ii = 1:length(max_clique)
    temp_clique = max_clique{ii}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Chordal relaxation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display('--------- Chordal calculation ----------')

cvx_begin
    variables Pg(n) Qg(n) Pinj(n) Qinj(n) Vsq(n) aux(n); 
    variables Pf(m) Pt(m);
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
    
        
        Pg <= PgMax;
        Pg >= PgMin;
        Qg <= QgMax;
        Qg >= QgMin;
        Vsq >= WMin;
        Vsq <= WMax;        
                
        Pf <= Fmax;
        Pt <= Fmax;
        
        for i = 1:length(max_clique)
            temp_clique = max_clique{i};
            temp_matrix = W(temp_clique, temp_clique);
            temp_matrix == hermitian_semidefinite( length(temp_clique) );
        end
         
cvx_end
toc
objective_value = sum(aux)*conditionObj

maxEigRatio = 0;
for i = 1:length(max_clique)
        temp_clique = max_clique{i};
        temp_matrix = W(temp_clique, temp_clique);
        eig_lst = eigs(temp_matrix);
        maxEigRatio = max(maxEigRatio, eig_lst(2) / eig_lst(1));
end

maxEigRatio



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Matpower's solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display('--------- MATPOWER optimization ----------')
% opt = mpoption('OPF_FLOW_LIM', 1);
% results = runopf(mpc, opt);
% 
% obj = results.f;
% display(strcat('Total cost = ', num2str(obj)));


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

% % % % Create list of neighbours for each node
% % % % neighbours = {};
% % % % for i = 1:n
% % % %     neighbours{i} = {};
% % % % end
% % % % 
% % % % for i = 1:length(mpc.branch(:, 1))
% % % %     len_1 = length(neighbours{mpc.branch(i, 1)});
% % % %     len_2 = length(neighbours{mpc.branch(i, 2)});
% % % %     neighbours{mpc.branch(i, 1)}{len_1 + 1} = mpc.branch(i, 2);
% % % %     neighbours{mpc.branch(i, 2)}{len_2 + 1} = mpc.branch(i, 1);
% % % % end
