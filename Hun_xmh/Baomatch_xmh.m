function [assign_mat, cost_tol] = Baomatch_xmh(cost_mat)
% This function is 
% Input:
% Output:

cost_matrix = cost_mat;

% get the number of columns and rows of the matrix
[rows, cols] = size(cost_matrix);
flag_t = 0;
if rows < cols
    cost_matrix = cost_matrix';
    [rows, cols] = size(cost_matrix);
    flag_t = 1;
end

cost_tol = inf;
assign_mat = zeros(rows, 2);
assign_mat(:, 1) = 1 : rows;
perms_all = perms(1 : rows);

for m_idx = 1 : size(perms_all, 1)
    match_s = perms_all(m_idx, :);
    match_s(match_s > cols) = 0;
    % 
    cost_s = 0;
    for r_idx = 1 : rows
        if match_s(r_idx) == 0
            continue;
        end
        cost_s = cost_s + cost_matrix(r_idx, match_s(r_idx));
    end
    if cost_s < cost_tol
        assign_mat(:, 2) = match_s';
        cost_tol = cost_s;
    end
end

assign_mat(assign_mat(:, 2) == 0, :) = [];
if flag_t == 1
    assign_mat = flip(assign_mat, 2);
end
assign_mat = sortrows(assign_mat, 1);
cost_tol = ...
    sum(cost_mat(sub2ind(size(cost_mat), assign_mat(:, 1), assign_mat(:, 2))));

end