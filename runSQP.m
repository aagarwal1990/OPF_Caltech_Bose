function [ V_fin ] = runSQP( V0, case_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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



% new code to determine sj
from_verts = mpc.branch(:,1);
to_verts = mpc.branch(:,2);

res = mpc.branch(:,3);
rctn = mpc.branch(:,4);


Adj = zeros(n);
Adm = zeros(n)
for ii = 1:m
    Adj(from_verts(ii), to_verts(ii)) = 1;
    Adj(to_verts(ii), from_verts(ii)) = 1;
    
    y_jk = complex(res(ii), -rctn(ii))/sqrt(res(ii), rctn(ii))
    Adm(from_verts(ii), to_verts(ii)) = y_jk;
    Adm(to_verts(ii), from_verts(ii)) = y_jk;
end

% to ensure that diagonal entries are zero
Adj(logical(eye(n)))=0;


% iterative step to determine s for each iteration
% assuming V_k has been determined

% matrix containing jth row multiplied by (V_k)_j
A = Adj * diag(V_k);
B = ((temp-temp').*Adm) * diag(V_k);
s = B(logical(eye(n)));

end

