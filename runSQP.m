function [ V_fin ] = runSQP( V0, case_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new code to determine sj
from_verts = mpc.branch(:,1);
to_verts = mpc.branch(:,2);

res = mpc.branch(:,3);
rctn = mpc.branch(:,4);


Adj = sparse(n,n);
Admittance = sparse(n,n)
for ii = 1:m
    jj = from_verts(ii);
    kk = to_verts(ii);
    Adj(jj, kk) = 1;
    Adj(kk, jj) = 1;
    
    y_jk = complex(res(ii), -rctn(ii))/sqrt(res(ii)^2 + rctn(ii)^2)
    Admittance(jj, kk) = y_jk;
    Admittance(kk, jj) = y_jk;
end

% to ensure that diagonal entries are zero
Adj(logical(eye(n)))=0;

rAdmittance = real(Admittance);% real part of admittance matrix
iAdmittance = imag(Admittance);% imaginary part of admittance matrix

nbrs = cell(1,n);
for jj = 1:n
    nbrs{jj} = find(Adj(jj,:));
end


% Following is code for iterative step which takes as input V_k
% iterative step to determine s for each iteration
% assuming V_k (column vector) has been determined



% % matrix containing jth row multiplied by (V_k)_j
% A = Adj * diag(V_k);
% B = ((temp-temp').*Admittance) * diag(V_k);
% s = B(logical(eye(n))); % this a part of g(x_k) in the book chapter 4
% 
% num_ineq = 2*(n+m); % min max for voltages and branch powers
% 
% lambda = zeros(num_ineq, 1);




% To determine the cost function which is dependent only on the real part
% as is seen in the solveOPF cvx, since cost is function of Pg only

s = zeros(n,1);

for jj = 1:n
    s(jj) = abs(V_k(jj))^2*sum(Admittance(:,jj)) - V_k(jj)*sum(Admittance(jj,:)*V_k);
end
Pinj = real(s);
Qinj = imag(s);

Vsq = abs(V_k).^2;
rV_k = real(V_k);% real part of voltage vector
iV_k = imag(V_k);% imaginary part of voltage vector
% no we have the constant parts of the inequalities . Need to evaluate the linear parts
% for that we use jacobian matrix in the various components

% To determine the jacobian matrices
jacobian_Pinj = sparse(2*n, n);
jacobian_Qinj = sparse(2*n, n);
for jj = 1:n
    
    jacobian_Pinj(2*jj-1,jj) = 2*rV_k(jj)*sum(rAdmittance(jj,:)) ...
        + rAdmittance(jj,:)*rV_k ...
        - iV_k(jj)*sum(iAdmittance(jj,:)) ...
        + iAdmittance(jj,:)*iV_k ...
        + iV_k(jj)*rAdmittance(JJ,:)*iV_k; % d /d real()
    jacobian_Pinj(2*jj, jj) = -2*iV_k(jj)*rAdmittance(jj,:)*iV_k ...
        + rV_k(jj)*iAdmittance(jj,:)*rV_k ...
        - rV_k(jj)*iAdmittance(jj:1)*i
    jacobian_Qinj
    jacobian_Qinj
    
    nbrs_jj = nbrs{jj}
    for ii = 1:size(nbrs_jj, 2)
        jacobian_Pinj(2*jj-1, ii) = 
        jacobian_Pinj(2*jj, ii) = 
        jacobian_Qinj
        jacobian_Qinj
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% this is to be done multipe times
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
        
%         % Line limits
%         for bb = 1:m
%             Pf(bb) == real(trace(Ff{bb} * W));
%             Pt(bb) == real(trace(Tt{bb} * W));
%         end
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

end

