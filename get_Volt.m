function [ Volt, thetas ] = get_Volt(from_list, to_list, num_verts, W)
% Function takes in W matrix and edge list. Outputs list of  voltages 
% for each node in network

% Volt is a list of complex values of voltages for each node. Magnitude of
% voltage of node i given by sqrt of magnitude of Wii. Angle given by 
% theta_i = theta_j - angle btw. Wij where theta_j is known
Volt = zeros(num_verts, 1);

% list of thetas initialized to -1. theta for node 1 set to 0
thetas = [];
for node = 1:num_verts
    thetas(node) = -1;
end
thetas(1) = 0;
% Add voltage to Volt vector
magn_1 = sqrt(abs(W(1, 1)));
Volt(1) = magn_1 * exp(1i * 0);

% Create an adjacency matrix of the spanning tree of power network
adjacency_matrix = sparse(get_adj_matrix(from_list, to_list, num_verts));
span_tree = graphminspantree(adjacency_matrix);
span_tree_adj_matrix = span_tree + span_tree';


% get order of traversing graph through a BFS of spanning tree adjMatrix
order = graphtraverse(span_tree_adj_matrix, 1, 'Method', 'BFS');
% list of nodes for which we do not know theta yet
unchecked_nodes = 2:num_verts;

for ii = 1:length(order)
    node = order(ii);
    % find all nodes adjacent to current node being considered in BFS
    adj_nodes = graphtraverse(span_tree_adj_matrix, node, 'depth', 1);
    % find intersection of unchecked nodes and adjacent nodes
    unchecked_child_nodes = intersect((adj_nodes), (unchecked_nodes));
    % iterate through unchecked adjacent nodes and find angle
    for jj = 1:length(unchecked_child_nodes)
        temp_node = unchecked_child_nodes(jj);
        % find relative angle between child node and current node 
        thetas(temp_node) = thetas(node) - angle(W(node, temp_node));
        % node has now been checked, remove it from unchecked nodes vector
        unchecked_nodes = unchecked_nodes(unchecked_nodes ~= temp_node);
        % Add voltage to Volt vector
        magn = sqrt(abs(W(temp_node, temp_node)));
        ang = thetas(temp_node);
        Volt(temp_node) = magn * exp(1i * ang);
    end
end


        