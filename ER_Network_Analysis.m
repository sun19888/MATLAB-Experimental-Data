% Set parameters
n1 = 50; % Number of nodes
n2 = 100;
n3 = 300;
p = 0:0.02:1;   % Edge probability
num_iterations = 50; % Number of iterations

% Initialize arrays to store all results
all_node1 = zeros(num_iterations, length(p));
all_edge1 = zeros(num_iterations, length(p));
all_node2 = zeros(num_iterations, length(p));
all_edge2 = zeros(num_iterations, length(p));
all_node3 = zeros(num_iterations, length(p));
all_edge3 = zeros(num_iterations, length(p));

% Calculate node and edge metrics
for i = 1:length(p)
    for k = 1:num_iterations
        [all_node1(k,i), all_edge1(k,i)] = ernet(n1, p(i));
        [all_node2(k,i), all_edge2(k,i)] = ernet(n2, p(i));
        [all_node3(k,i), all_edge3(k,i)] = ernet(n3, p(i));
    end
end

% Calculate mean and standard deviation
node_values1 = mean(all_node1);
edge_values1 = mean(all_edge1);
node_std1 = std(all_node1);
edge_std1 = std(all_edge1);

node_values2 = mean(all_node2);
edge_values2 = mean(all_edge2);
node_std2 = std(all_node2);
edge_std2 = std(all_edge2);

node_values3 = mean(all_node3);
edge_values3 = mean(all_edge3);
node_std3 = std(all_node3);
edge_std3 = std(all_edge3);

% Create x values
x = p;

% Plot curves
figure;
hold on;

% Plot node metric curves with error bars
errorbar(x, node_values1, node_std1, 'k-', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x, node_values2, node_std2, 'r--', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x, node_values3, node_std3, 'b-.', 'LineWidth', 1.2, 'CapSize', 5);

% Plot edge metric curves with error bars
errorbar(x, edge_values1, edge_std1, 'k-s', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x, edge_values2, edge_std2, 'r--o', 'LineWidth', 1.2, 'CapSize', 5);
errorbar(x, edge_values3, edge_std3, 'b-.*', 'LineWidth', 1.2, 'CapSize', 5);

% Set Y-axis display range
ylim([-0.02, 1]);

% Add title and labels
xlabel('Edge Probability (p)');
ylabel('Node/Edge Metric (R_n/R_e)');
legend('R_n_5_0', 'R_n_1_0_0', 'R_n_3_0_0', 'R_e_5_0', 'R_e_1_0_0', 'R_e_3_0_0', ...
       'FontSize', 12, 'Location', 'best');

hold off;

% Calculate ER network metrics
function [avg_node, avg_edge] = ernet(n, p)
    % Initialize cumulative values
    sum_node = 0;
    sum_edge = 0;
    
    % Loop 50 times
    num_iterations = 50;
    for k = 1:num_iterations
        A = generate_ER_graph(n, p);
        numEdges = round((nnz(A) - nnz(diag(A))) / 2); % Count undirected edges
        
        % Calculate structural controllability metrics
        U = 0;
        C = 0;
        D = size(A);
        E = eig(A);
        for i = 1:D(1)
            C = D(1) - rank(E(i) * eye(D(1)) - A);
            if C >= U
                U = C;
            end
        end
        
        % Calculate node and edge metrics and accumulate
        node = U / n;
        edge = ((n - U) + (U - 1)) / numEdges;
        if edge > 1
            edge = 1; % Cap edge metric at 1
        end
        sum_node = sum_node + node;
        sum_edge = sum_edge + edge;
    end
    
    % Calculate averages
    avg_node = sum_node / num_iterations;
    avg_edge = sum_edge / num_iterations;
end

% Generate ER random graph adjacency matrix
function adjacency_matrix = generate_ER_graph(n, p)
    adjacency_matrix = zeros(n);
    for i = 1:n
        for j = i+1:n
            if rand() < p
                adjacency_matrix(i, j) = 1;
                adjacency_matrix(j, i) = 1; % Undirected edge
            end
        end
    end
end