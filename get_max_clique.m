function [ max_cliques ] = get_max_clique( nbrs, sigma, inv_sigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = length(nbrs);

% for all verts of G do
size_v = zeros(1,n);
max_cliques = {};
num_max_cliques = 0;

% for all i from 1 to n dp
for ii = 1:n
    vert = inv_sigma(ii);
    if length(nbrs{vert}) == 0
        num_max_cliques = num_max_cliques + 1;
        max_cliques{num_max_cliques} = [vert];
    end
    
    % set X to {w belong to..}
    X = {};
    X_count = 0;
    nbrs_v_list = nbrs{vert};
    for jj = 1:length(nbrs_v_list)
        nbr_v = nbrs_v_list(jj);
        if sigma(vert) < sigma(nbr_v)
            X_count = X_count + 1;
            X{X_count} = nbr_v;
        end
    end
    
    X = cell2mat(X);

    if length(X) ~= 0
        if size_v(vert) < length(X)
            num_max_cliques = num_max_cliques + 1;
            max_cliques{num_max_cliques} = cat(2,[vert],X);
        end
        
        m_v = inv_sigma(min(sigma(X)));
        size_v(m_v) = max(size_v(m_v), length(X)-1);
    end
end

display('max_cliques_concise');
max_len_clique = 0;
for ii = 1:length(max_cliques)
    if max_len_clique < length(max_cliques{ii})
        max_len_clique = length(max_cliques{ii});
    end
end
M = zeros(num_max_cliques, max_len_clique);
for ii = 1:length(max_cliques)
    M(ii, 1:length(max_cliques{ii})) = max_cliques{ii};
end

