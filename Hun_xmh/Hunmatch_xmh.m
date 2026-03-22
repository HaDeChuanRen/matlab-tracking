function [assign_mat, cost_tol] = Hunmatch_xmh(cost_mat)
% This function is 
% Input:
% Output:

cost_matrix = cost_mat;
cost_tol = 0;
% get the number of columns and rows of the matrix
[rows, cols] = size(cost_matrix);
% cost_matrix = cost_matrix - min(cost_matrix(:));

% 
% if rows >= cols
%     cost_matrix = [cost_matrix, inf(rows, rows - cols)];
%     min_rc = cols;
%     cols = rows;
% elseif rows < cols
%     cost_matrix = [cost_matrix; inf(cols - rows, cols)];
%     min_rc = rows;
%     rows = cols;
% end

% initialization
assign_mat = zeros(min(rows, cols), 2);
% assign_mat(:, 1) = 1 : rows;
% marked_rows = false(rows, 1);
% marked_cols = false(cols, 1);


% rows minus
% for r_idx = 1 : rows
%     min_ral = min(cost_matrix(r_idx, :));
%     cost_matrix(r_idx, :) = cost_matrix(r_idx, :) - min_ral;
% end

% columns minus
% for c_idx = 1 : cols
%     min_cal = min(cost_matrix(:, c_idx));
%     cost_matrix(:, c_idx) = cost_matrix(:, c_idx) - min_cal;
% end

% initialization
M_match = zeros(0, 2);
y_pot = zeros(rows + cols, 1);
Gmat_y = false(rows, cols);
while 1
    Rest_r = setdiff(1 : rows, M_match(:, 1));
    Rest_c = setdiff(rows + 1 : rows + cols, M_match(:, 2));
    
    if isempty(Rest_r) || isempty(Rest_c)
        break;
    end
    Ymat_pot = y_pot(1 : rows) + y_pot(rows + 1 : end)';
    Gmat_y = abs(cost_matrix - Ymat_pot) < 1e-6;
%     Gmat_y(abs(cost_matrix - Ymat_pot) < 1e-10) = 1;
%     Gmat_y(abs(cost_matrix - Ymat_pot) > 1e-10) = 0;
    [Z_ver, P_paths] = DFS_Hun(Rest_r, Rest_c, Gmat_y, M_match);
    if ~isempty(P_paths)
        c_insc = intersect(Z_ver, Rest_c);
        for p_idx = 1 : length(P_paths)
            if P_paths{p_idx}(end) == c_insc(1)
                p_path = P_paths{p_idx};
                for n_idx = 1 : 2 : (length(p_path) - 1)
                    M_match(M_match(:, 2) == (p_path(n_idx + 1)), :) = [];
                    M_match(end + 1, :) = [p_path(n_idx),...
                        p_path(n_idx + 1)];
                end
                break;
            end
        end
    else
        r_set = intersect(Z_ver, 1 : rows);
        c_set = setdiff(rows + 1 : rows + cols, Z_ver);
        delta_all = zeros(length(r_set), length(c_set));
        for rs_idx = 1 : length(r_set)
            for cs_idx = 1 : length(c_set)
                r_ele = r_set(rs_idx);
                c_ele = c_set(cs_idx);
                delta_all(rs_idx, cs_idx) = cost_matrix(r_ele, ...
                    c_ele - rows) - y_pot(r_ele) - y_pot(c_ele);
            end
        end
        delta_now = min(delta_all(:));
        cs_set = intersect(rows + 1 : rows + cols, Z_ver);
        y_pot(r_set) = y_pot(r_set) + delta_now;
        y_pot(cs_set) = y_pot(cs_set) - delta_now;
        y_pos = sum(y_pot);
    end
end

% order M_match to form assign_mat
% for r_idx = 1 : rows
%     assign_mat(r_idx, 2) = M_match(M_match(:, 1) == r_idx, 2) - rows;
%     cost_tol = cost_tol + cost_mat(r_idx, assign_mat(r_idx, 2));
% end
M_match(:, 2) = M_match(:, 2) - rows;
assign_mat = sortrows(M_match, 1);
cost_tol = ...
    sum(cost_mat(sub2ind(size(cost_mat), assign_mat(:, 1), assign_mat(:, 2))));

end