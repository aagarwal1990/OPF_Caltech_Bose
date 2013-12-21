function [ M ] = get_adj_matrix( from_list, to_list, num_verts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = num_verts;
m = length(from_list);

M = zeros(n);

for i = 1:m
    M(from_list(i), to_list(i)) = 1;
    M(to_list(i), from_list(i)) = 1;
end

end

