function [ cell_c, M_ch ] = get_chordal_matrix_chole_random( from_list, to_list, n, num_iter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

M_adj = get_adj_matrix(from_list, to_list, n);

D = diag(sum(M_adj));
L = sparse(D - M_adj);

min_fill = n^2;
min_fill_list = zeros(1, num_iter);

for ii = 1:num_iter
    
    P = sparse(zeros(n));
    permute_list = randperm(n);
    for jj = 1:n
        P(jj, permute_list(jj)) = 1;
    end

    temp = P*(L+eye(size(L)))*P';
    R = chol(temp);

    M = sparse(triu(R,1)'+R);
    nnz_M = nnz(M);
    if nnz_M < min_fill
        min_fill = nnz_M;
        min_fill_list(ii) = min_fill;
        M_ch = P'*M*P;
    end
end    
cell_c = {};
for i = 1:n
    temp = find(M_ch(:,i))';
    cell_c{i} = temp(temp ~= i);
end

end
