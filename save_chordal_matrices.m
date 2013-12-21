case_list = [9,14,30,39,57];


for ii = 1:length(case_list)
    case_num = case_list(ii);
    case_name = strcat('case', int2str(case_num));
    mpc = loadcase(case_name);

    from_list = mpc.branch(:,1);
    to_list = mpc.branch(:,2);
    n = length(mpc.bus(:,1));
    m = length(mpc.branch(:,1));
    
    [nbrs, M_ch] = get_chordal_matrix_chole_random(from_list, to_list, n, ceil(n*n*10));
    
    [sigma, PEO] = get_PEO(nbrs);

    % Get Max Cliques
    max_clique = get_max_clique(nbrs, sigma, PEO);
    
    file_name = strcat('chordal_matrices/ch_mat_', case_name, '.mat');
    save(file_name, 'M_ch','nbrs', 'max_clique');
end