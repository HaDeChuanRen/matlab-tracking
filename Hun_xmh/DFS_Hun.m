function [Z_ver, P_paths] = DFS_Hun(Rest_r, Rest_c, Gmat_y, M_match)
% breadth-first search

P_paths = {};
Z_ver = zeros(0, 1);
rows = size(Gmat_y, 1);
cols = size(Gmat_y, 2);
Num_ver = rows + cols;
GM_graph = false(Num_ver);

G_all = find(Gmat_y == 1);

for g_idx = 1 : length(G_all)
    [rnode_g, cnode_g] = ind2sub([rows, cols], G_all(g_idx));
    GM_graph(rnode_g, cnode_g + rows) = 1;
end

for m_idx = 1 : size(M_match, 1)
    rnode_m = M_match(m_idx, 1);
    cnode_m = M_match(m_idx, 2);
    GM_graph(rnode_m, cnode_m) = 0;
    GM_graph(cnode_m, rnode_m) = 1;
end

for r_idx = Rest_r
    visited_node = false(Num_ver, 1);
    visited_node(r_idx) = 1;
    queue_r = r_idx;
    flag_depth = 1;
    queue_depth = flag_depth;
    p_path = [];
    while ~isempty(queue_r)
        current_node = queue_r(1);
        p_path = [p_path, current_node];
        Z_ver = unique([Z_ver; current_node]);
        queue_r(1) = [];
        queue_depth(1) = [];
        neigh_all = find(GM_graph(current_node, :));
        neigh_all(visited_node(neigh_all)) = [];
        if isempty(neigh_all)
            if sum(ismember(Rest_c, current_node))
                P_paths{end + 1} = p_path;
                % Z_ver = setdiff(Z_ver, Rest_r);
                % return;
            end
            if ~isempty(queue_depth)
                delta_depth = flag_depth - queue_depth(1);
                p_path(end - delta_depth : end) = [];
                flag_depth = queue_depth(1);
            else
                p_path = [];
                flag_depth = 1;
            end
        else
            flag_depth = flag_depth + 1;
        end
        for neigh_node = neigh_all
            visited_node(neigh_node) = 1;
            queue_r = [neigh_node, queue_r];
            
            queue_depth = [flag_depth, queue_depth];
        end
    end
end
% Z_ver = setdiff(Z_ver, Rest_r);

end



% if current_node <= rows
%     neigh_node(visited_node(neigh_node)) = [];
% else
%     neigh_node = M_match(M_match(:, 2) == (current_node - rows));
%     neigh_node(visited_node(neigh_node)) = [];
% end

% for r_idx = Rest_r
%     visited_node = false(Num_ver, 1);
%     visited_node(r_idx) = 1;
%     queue_r = r_idx;
%     %¡¡flag = 1;
%     while ~isempty(queue_r)
%         current_node = queue_r(1);
%         queue_r(1) = [];
%         p_path = current_node;
%         queue_c = find(Gmat_y(current_node, :)) + rows;
%         queue_c(visit(queue_c)) == [];
% 
%         while 1
%             
%             if ~isempty(neigh_node)
%                 visited_node(neigh_node) = 1;
%                 queue_r = [queue_r, neigh_node];
%                 p_path = [r_idx, queue_r];
%                 Z_ver = unique([Z_ver; neigh_node]);
%             else
%                 P_paths{end + 1} = p_path;
%                 break;
%             end
%         end
%     end
% end