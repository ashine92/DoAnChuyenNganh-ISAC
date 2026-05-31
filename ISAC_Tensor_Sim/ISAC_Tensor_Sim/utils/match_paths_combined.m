function perm = match_paths_combined(az_est, el_est, tau_est, az_true, el_true, tau_true)
% Match estimated paths to true paths using combined cost on angles and ToA.
% 
% FIX (ITERATION 3 Phase 2b): Uses more robust matching algorithm:
% - Tries multiple matching strategies
% - Uses Hungarian algorithm for global optimality
% - Falls back to greedy if needed
% - Implemented to fix path reordering issue
%
% Input:
%   az_est, el_est, tau_est   - Estimated angles and ToA (L×1)
%   az_true, el_true, tau_true - Ground truth (L×1)
%
% Output:
%   perm - Permutation vector (L×1), perm(i) = j means est_i matches true_j

    L = length(az_true);
    if L == 1
        perm = 1;
        return;
    end

    % Normalized cost matrix (combine angle and ToA errors)
    cost_az  = abs(az_est(:) - az_true(:)');
    cost_el  = abs(el_est(:) - el_true(:)');
    cost_tau = abs(tau_est(:) - tau_true(:)');

    % Normalize each cost by its range for balanced weighting
    range_az  = max(eps, max(cost_az(:)));
    range_el  = max(eps, max(cost_el(:)));
    range_tau = max(eps, max(cost_tau(:)));

    % Combined normalized cost
    cost_combined = cost_az/range_az + cost_el/range_el + cost_tau/range_tau;

    % ===== ROBUST MATCHING: Try multiple strategies =====
    
    % Strategy 1: Greedy matching (original)
    perm_greedy = greedy_match(cost_combined);
    
    % Strategy 2: Hungarian algorithm (optimal for L×L assignment)
    % Try to use Hungarian if available, otherwise greedy
    perm_hungarian = [];
    try
        % Check if optimization toolbox is available
        assignmentproblem = @(c) hungarianMatch(c);
        perm_hungarian = hungarianMatch(cost_combined);
    catch
        % Hungarian not available, will use greedy
    end
    
    % Strategy 3: Closest-by-distance greedy (less ambiguity)
    perm_closest = closest_match(cost_combined);
    
    % ===== SELECT BEST PERMUTATION =====
    % Compute "confidence" for each permutation (lower cost = higher confidence)
    cost_greedy = trace_cost(cost_combined, perm_greedy);
    cost_closest = trace_cost(cost_combined, perm_closest);
    
    if ~isempty(perm_hungarian)
        cost_hungarian = trace_cost(cost_combined, perm_hungarian);
        
        % Use Hungarian if it's better
        if cost_hungarian < min(cost_greedy, cost_closest)
            perm = perm_hungarian;
            return;
        end
    end
    
    % Otherwise use the better of greedy vs closest
    if cost_closest < cost_greedy
        perm = perm_closest;
    else
        perm = perm_greedy;
    end
end


function perm = greedy_match(cost_mat)
    % Original greedy matching
    L = size(cost_mat, 1);
    perm = zeros(1, L);
    used = false(1, L);
    
    for i = 1:L
        c = cost_mat(i,:);
        c(used) = inf;
        [~, best_j] = min(c);
        perm(i) = best_j;
        used(best_j) = true;
    end
end


function perm = closest_match(cost_mat)
    % Match paths by closest distance (min over all)
    % This prioritizes globally closest matches
    L = size(cost_mat, 1);
    perm = zeros(1, L);
    used_est = false(1, L);
    used_true = false(1, L);
    
    % Create list of all costs with indices
    cost_list = [];
    for i = 1:L
        for j = 1:L
            cost_list = [cost_list; cost_mat(i, j), i, j];
        end
    end
    
    % Sort by cost (ascending)
    [~, sort_idx] = sort(cost_list(:, 1));
    cost_list = cost_list(sort_idx, :);
    
    % Greedily match lowest-cost pairs
    matches = 0;
    for idx = 1:size(cost_list, 1)
        if matches == L
            break;
        end
        
        cost_val = cost_list(idx, 1);
        i_est = cost_list(idx, 2);
        j_true = cost_list(idx, 3);
        
        if ~used_est(i_est) && ~used_true(j_true)
            perm(i_est) = j_true;
            used_est(i_est) = true;
            used_true(j_true) = true;
            matches = matches + 1;
        end
    end
end


function perm = hungarianMatch(cost_mat)
    % Hungarian algorithm for optimal assignment
    % Returns permutation that minimizes total cost
    %
    % Simple recursive implementation for small L (≤ 3)
    
    L = size(cost_mat, 1);
    
    if L == 1
        perm = 1;
        return;
    end
    
    if L == 2
        % Try both permutations
        cost1 = cost_mat(1,1) + cost_mat(2,2);
        cost2 = cost_mat(1,2) + cost_mat(2,1);
        
        if cost1 <= cost2
            perm = [1, 2];
        else
            perm = [2, 1];
        end
        return;
    end
    
    if L == 3
        % Try all 6 permutations for L=3
        perms = [1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1];
        costs = zeros(6, 1);
        
        for p = 1:6
            cost_p = 0;
            for i = 1:3
                cost_p = cost_p + cost_mat(i, perms(p, i));
            end
            costs(p) = cost_p;
        end
        
        [~, best_p] = min(costs);
        perm = perms(best_p, :);
        return;
    end
    
    % For larger L, fall back to greedy
    perm = greedy_match(cost_mat);
end


function total_cost = trace_cost(cost_mat, perm)
    % Compute total cost for a given permutation
    L = length(perm);
    total_cost = 0;
    for i = 1:L
        total_cost = total_cost + cost_mat(i, perm(i));
    end
end

