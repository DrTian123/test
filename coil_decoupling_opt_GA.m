%% MATLAB Code: Optimize Coil Spacing for Minimum Mutual Inductance
% Two adjacent rectangular planar coils

% Parameters
coil_width = 0.4;      % Coil width in meters (40 cm)
coil_height = 0.6;     % Coil height in meters (60 cm)
num_turns = 10;        % Number of turns per coil
distance = 0.1;        % Distance between the two coils (meters)

% Genetic Algorithm Parameters
population_size = 50;  % Number of individuals in the population
max_generations = 100; % Maximum number of generations

% Define the objective function
% The goal is to minimize the mutual inductance between the two coils
obj_fun = @(x) mutual_inductance(x, coil_width, coil_height, distance, num_turns);

% Define bounds for optimization
lb = zeros(1, num_turns - 1); % Lower bounds (minimum spacing between turns = 0)
ub = ones(1, num_turns - 1) * (coil_height / (num_turns - 1)); % Upper bounds (maximum spacing)

% Run the genetic algorithm
options = optimoptions('ga', 'PopulationSize', population_size, ...
                        'MaxGenerations', max_generations, ...
                        'Display', 'iter');

[optimal_spacing, min_mutual_inductance] = ga(obj_fun, num_turns - 1, [], [], [], [], lb, ub, [], options);

% Display results
disp('Optimal spacing between turns:');
disp(optimal_spacing);
disp('Minimum mutual inductance:');
disp(min_mutual_inductance);

%% Objective Function: Mutual Inductance Calculation
function M = mutual_inductance(spacing, width, height, dist, turns)
    % Add fixed spacing for the first turn
    spacing = [0, cumsum(spacing)];

    % Ensure total height constraint
    if spacing(end) > height
        M = 1e6; % Penalize invalid configurations
        return;
    end

    % Define coil positions for each turn
    y_positions = spacing;

    % Calculate mutual inductance between two coils
    M = 0;
    for i = 1:turns
        for j = 1:turns
            % Approximate mutual inductance between two rectangular loops
            M = M + rectangular_loop_mutual(width, y_positions(i), y_positions(j), dist);
        end
    end
end

%% Function: Mutual Inductance Between Two Rectangular Loops
function M = rectangular_loop_mutual(width, y1, y2, distance)
    % Approximate the mutual inductance between two parallel rectangular loops
    % Using the Neumann formula approximation
    mu0 = 4 * pi * 1e-7; % Permeability of free space

    % Distance between the centers of two loops
    d = sqrt((y1 - y2)^2 + distance^2);

    % Simplified mutual inductance formula (approximation for distant loops)
    M = mu0 * width^2 / (2 * pi * d);
end
