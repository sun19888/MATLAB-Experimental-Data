% Set parameters
n1 = 50; % Number of nodes
n2 = 100;
n3 = 300;
m = 0:1:10; % Number of new edges added each time a new node is added
num_iterations = 50; % Define number of iterations

% Initialize arrays to store node and edge metrics
node_values1 = zeros(size(m));
edge_values1 = zeros(size(m));
node_values2 = zeros(size(m));
edge_values2 = zeros(size(m));
node_values3 = zeros(size(m));
edge_values3 = zeros(size(m));

% Initialize arrays to store all iteration results
all_node1 = zeros(num_iterations, length(m));
all_edge1 = zeros(num_iterations, length(m));
all_node2 = zeros(num_iterations, length(m));
all_edge2 = zeros(num_iterations, length(m));
all_node3 = zeros(num_iterations, length(m));
all_edge3 = zeros(num_iterations, length(m));

% Calculate node and edge metrics
for i = 1:length(m)
    for iter = 1:num_iterations
        [all_node1(iter, i), all_edge1(iter, i)] = banet(n1, m(i));
        [all_node2(iter, i), all_edge2(iter, i)] = banet(n2, m(i));
        [all_node3(iter, i), all_edge3(iter, i)] = banet(n3, m(i));
    end
    
    node_values1(i) = mean(all_node1(:, i));
    edge_values1(i) = mean(all_edge1(:, i));
    node_values2(i) = mean(all_node2(:, i));
    edge_values2(i) = mean(all_edge2(:, i));
    node_values3(i) = mean(all_node3(:, i));
    edge_values3(i) = mean(all_edge3(:, i));
end

% Calculate standard deviation
node_std1 = std(all_node1);
edge_std1 = std(all_edge1);
node_std2 = std(all_node2);
edge_std2 = std(all_edge2);
node_std3 = std(all_node3);
edge_std3 = std(all_edge3);

% Create x values
x = m;

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

% Set display range for X and Y axes
ylim([-0.02, 1]);

% Add title and labels
xlabel('<k>/2');
ylabel('R_n/R_e');
legend('R_n_5_0', 'R_n_1_0_0', 'R_n_3_0_0', 'R_e_5_0', 'R_e_1_0_0', 'R_e_3_0_0', ...
       'FontSize', 12, 'Location', 'best');
hold off;

% Calculate BA network metrics
function [avg_node, avg_edge] = banet(n, m)
    % Create adjacency matrix
    A = zeros(n);
    
    % Initialize network with m+1 nodes
    A(1:m+1, 1:m+1) = randi([0, 1], m+1, m+1);
    
    % Initialize cumulative variables
    total_node = 0;
    total_edge = 0;
    
    % Loop 50 times
    num_iterations = 50;
    
    for iter = 1:num_iterations
        % Add edges for new nodes
        for newNode = m+2:n
            % Calculate degrees of existing nodes
            degrees = sum(A, 2);
            
            % Calculate connection probabilities
            probabilities = degrees / sum(degrees);
            
            % Randomly select m nodes from probability distribution
            selectedNodes = [];
            for i = 1:m
                selectedNode = randsample(newNode-1, 1, true); 
                selectedNodes = [selectedNodes, selectedNode];
            end
            
            % Add edges
            A(newNode, selectedNodes) = 1;
            A(selectedNodes, newNode) = 1;
        end
        
        % Calculate maximum geometric multiplicity
        U = 0;
        E = eig(A);
        
        for i = 1:length(E)
            C = n - rank(E(i) * eye(n) - A);
            if C >= U
                U = C;
            end
        end
        
        % Calculate number of nodes and edges
        numNodes = size(A, 1);
        numEdges = round((nnz(A) - nnz(diag(A))) / 2); % Round edge count
        
        % Calculate node and edge metrics
        node = U / numNodes;
        edge = ((n - U) + (U - 1)) / numEdges;
        
        if edge > 1
            edge = 1;
        end
        
        % Accumulate node and edge metric values
        total_node = total_node + node;
        total_edge = total_edge + edge;   
        
        % Reset adjacency matrix
        A = zeros(n);
        A(1:m+1, 1:m+1) = 1;
    end 
    
    % Calculate averages
    avg_node = total_node / num_iterations;
    avg_edge = total_edge / num_iterations;
end