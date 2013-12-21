function [ sigma, inv_sigma ] = get_PEO( nbrs )
% INPUT -

% import java.util.LinkedList
% 
% list = java.util.LinkedList;
% PEO_list = java.util.LinkedList;
% n = length(M_ch);
% label = n;
% 
% % select random vertex as curr vertex
% curr_vert = randi(n);
% PEO_list.add(curr_vert);
% % update neighbours and add to unvisited list
% nbrs = find(M_ch(:,curr_vert));
% for ii = 1:length(nbrs)
%     list.add({nbrs(ii), [label], 1});
% end
% 
% % select next vertex
% deci_power = @x 10.^x
% next_vert_index = -1;
% next_label_val;
% next_label_size = 0;
% for ii = 0:list.size()-1
%     temp = cell(list.get(ii));
%     temp_vert = temp{1};
%     temp_label = temp{2};
%     temp_label_size; = temp{3};
%     temp_label_val = sum(temp_label.*deci_power(flipr[1:temp_n_label]))
%     if next_label_val < temp_label_val
%         next_vert_index = ii;
%     end
%     next_vert = 
% end
% % make next vertex current vertex
% % update neighbours
% % delete curr vertex from list of unvisited vertices


n = length(nbrs);
% for all verts set label = 0
label_list = zeros(1,n);

sigma = zeros(1,n);
inv_sigma = zeros(1,n);

for ii = n:-1:1
    % choose an un numbered vertex v with largest label
    % works because all numbered vertices have label -1
    [val, vert] = max(label_list);
    sigma(vert) = ii;
    inv_sigma(ii) = vert;
    % making label for the vertex -1 because it has been numbered
    nbrs_v_list = nbrs{vert};
    % for all unnumbered vertices adjacent to vertex v
    for jj = 1:length(nbrs_v_list)
        w = nbrs_v_list(jj);
        if label_list(w) >= 0
            label_list(w) = label_list(w) + 1;
        end
    end
    label_list(vert) = -1;
end

end

