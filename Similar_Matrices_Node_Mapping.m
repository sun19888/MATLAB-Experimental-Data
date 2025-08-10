% System parameters
A3 = [
    0.56, -0.08, -0.16, -0.16,  1.04, -0.48,  0.08,  0.32;
   -6.08,  1.04,  2.88,  0.88, -0.72,  0.64, -0.44, -0.76;
   -5.84,  1.12,  3.24,  0.24,  1.44, -1.28, -0.12, -0.48;
   -2.48,  0.64,  1.28,  2.28, -0.32, -0.16, -0.64, -0.56;
   -1.12,  0.16,  0.32,  0.32, -0.08,  0.96, -0.16, -0.64;
   -8.92,  1.56,  4.12,  1.12, -0.28,  1.36, -0.56, -2.24;
    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00;
  -16.40,  3.20,  8.40,  2.40, -1.60, -0.80, -1.20, -2.80 
];
A3_prime = [
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 3, 0, 0, 3;
    0, 0, 1, 1, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 3, 0, 0;
    0, 0, 0, 0, 2, 0, 3, 0;
    0, 0, 0, 0, 0, 0, 0, 3;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 2  
];

[matchPairs, maxMatchNum] = similar_matrices_node_mapping(A3,A3_prime)

function [matchPairs, maxMatchNum] = similar_matrices_node_mapping(A, B)
% SIMILAR_MATRICES_NODE_MAPPING Finds node correspondence between similar matrices
% Input:
%   A, B - Adjacency matrices before and after similarity transformation
% Output:
%   matchPairs - Matching pairs [u v] where node u in transformed network corresponds to node v in original network
%   maxMatchNum - Number of maximum matched edges (checks if all nodes are matched)
    [P, R] = find_similarity_matrices(A, B);
    [matchPairs, maxMatchNum] = maxMatchingDirected(R)
end

function [P, R] = find_similarity_matrices(A, B)
% FIND_SIMILARITY_MATRICES Finds similarity transformation matrix and node correspondence matrix
% Input:
%   A - Adjacency matrix (n x n), directed graph (0-1 matrix)
% Output:
%   P - Similarity transformation matrix
%   R - Node correspondence matrix
    % Verify input matrices are square and same size
    if ~isequal(size(A), size(B)) || size(A,1) ~= size(A,2)
        error('Input matrices must be square matrices of same size');
    end

    [VA, DA] = eig(A);
    [VB, DB] = eig(B);
    
    P = VA / VB;
    
    P_inv = inv(P);
    
    n = size(A, 1);
    R = ones(n, n);
    
    % Populate R matrix
    for i = 1:n
        for j = 1:n
            if P_inv(i,j) * P(j,i) == 0
                R(i,j) = 0;
            end
        end
    end
end


function [matchPairs, maxMatchNum] = maxMatchingDirected(A)
% MAXMATCHINGDIRECTED Computes maximum matching for directed graph (structural controllability)
% Input:
%   A - Adjacency matrix (n x n), directed graph (0-1 matrix)
% Output:
%   matchPairs - Matching pairs [u v] indicating edge u->v is matched
%   maxMatchNum - Number of edges in maximum matching

    n = size(A,1);

    % 1. Construct bipartite graph (left=source nodes, right=target nodes)
    L = repmat((1:n)', 1, n);
    R = repmat(1:n, n, 1);
    edges = [L(A~=0), R(A~=0)]; % Bipartite edges: [left_node, right_node]

    % 2. Apply Hopcroft-Karp algorithm
    [pairU, pairV] = hopcroftKarp(n, n, edges);

    % 3. Format output
    matchPairs = [];
    for u = 1:n
        if pairU(u) ~= 0
            matchPairs = [matchPairs; u, pairU(u)];
        end
    end
    maxMatchNum = size(matchPairs, 1);
end

function [pairU, pairV] = hopcroftKarp(nU, nV, edges)
% HOPCROFTKARP Maximum bipartite matching using Hopcroft-Karp algorithm
% Input:
%   nU - Number of left vertices
%   nV - Number of right vertices
%   edges - E x 2 matrix where each row [u v] represents left u -> right v edge
% Output:
%   pairU - Right vertex matched to each left vertex
%   pairV - Left vertex matched to each right vertex

    % Build adjacency list
    adj = cell(nU,1);
    for i = 1:size(edges,1)
        u = edges(i,1);
        v = edges(i,2);
        adj{u} = [adj{u}, v];
    end
    pairU = zeros(nU,1);
    pairV = zeros(nV,1);
    dist  = zeros(nU,1);
    while bfs()
        for u = 1:nU
            if pairU(u) == 0
                dfs(u);
            end
        end
    end

    function ok = bfs()
        queue = [];
        for u = 1:nU
            if pairU(u) == 0
                dist(u) = 0;
                queue(end+1) = u;
            else
                dist(u) = inf;
            end
        end
        ok = false;
        qi = 1;
        while qi <= numel(queue)
            u = queue(qi); qi = qi + 1;
            for v = adj{u}
                if pairV(v) == 0
                    ok = true;
                elseif dist(pairV(v)) == inf
                    dist(pairV(v)) = dist(u) + 1;
                    queue(end+1) = pairV(v);
                end
            end
        end
    end

    function found = dfs(u)
        found = false;
        for v = adj{u}
            if pairV(v) == 0 || (dist(pairV(v)) == dist(u) + 1 && dfs(pairV(v)))
                pairU(u) = v;
                pairV(v) = u;
                found = true;
                return
            end
        end
        dist(u) = inf;
    end
end