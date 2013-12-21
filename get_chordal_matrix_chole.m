function [ cell_c, M ] = get_chordal_matrix_chole( from_list, to_list, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

M_adj = get_adj_matrix(from_list, to_list, n);

% P = zeros(n);
permute_list = 1:n;
% permute_list(2)=1;
% permute_list(1)=2;
permute_list(3)=4;
permute_list(4)=3;

for jj = 1:n
    P(permute_list(jj), jj) = 1;
end
sparse(P)

% P(1,1)=0;
% P(2,2)=0;
% P(1,2)=1;
% P(2,1)=1;

% p = symamd(M_adj);

M_adj = P*M_adj*P';
D = diag(sum(M_adj));
L = D - M_adj;
temp = L+eye(size(L));
R = chol(temp);
% R = chol(L+eye(size(L)));

M = triu(R,1)'+R;

M = P'*M*P;

cell_c = {};
for i = 1:n
    temp = find(M(:,i))';
    cell_c{i} = temp(temp ~= i);
end

end

