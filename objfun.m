function [objective, grad_cost] = objfun(exp_V_k)

case_num = 'case14';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optimization variables in real form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

for kk = 1:n
    Pinj(kk) = exp_V_k'* exp_Phi{kk} * exp_V_k ;
    Qinj(kk) = exp_V_k'* exp_Psi{kk} * exp_V_k;
    Vsq(kk)  = (exp_V_k(kk))^2 + (exp_V_k(kk + n))^2;
end

Pg = Pinj' + Pd;
Qg = Qinj' + Qd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write objective function and gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cost      = zeros(2*n,1);
grad_cost = zeros(2*n,1);
for kk = 1:n
    cost = cost + costGen1(kk)* Pg(kk);
    grad_cost = grad_cost + 2*costGen1(kk)*(exp_Phi{kk}*exp_V_k);
end
objective = sum(cost);

