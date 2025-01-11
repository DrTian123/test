% Optimization of coil spacing for rectangular planar coils decoupling
% This script optimizes the spacing between turns of adjacent rectangular coils
% to minimize mutual coupling

function [optimal_spacing] = coil_spacing_optimization()
    % Initial parameters
    coil_length = 0.1;  % Length of rectangular coil (m)
    coil_width = 0.1;   % Width of rectangular coil (m)
    wire_radius = 0.001; % Wire radius (m)
    num_turns = 5;      % Number of turns in each coil
    coil_separation = 0.02; % Vertical separation between coils (m)
    
    % Initial guess for spacing between turns
    initial_spacing = 0.005; % 5mm initial spacing
    
    % Optimization options
    options = optimset('Display', 'iter', 'TolX', 1e-6);
    
    % Perform optimization
    optimal_spacing = fminbnd(@(x) objective_function(x, coil_length, coil_width, ...
        wire_radius, num_turns, coil_separation), 0.001, 0.02, options);
    
    % Display results
    fprintf('Optimal spacing between turns: %.3f mm\n', optimal_spacing * 1000);
    
    % Plot the coils with optimal spacing
    plot_coils(optimal_spacing, coil_length, coil_width, num_turns);
end

function M = objective_function(spacing, length, width, wire_radius, num_turns, separation)
    % Calculate mutual inductance between two coils
    M = calculate_mutual_inductance(spacing, length, width, wire_radius, ...
        num_turns, separation);
    
    % Objective is to minimize absolute value of mutual inductance
    M = abs(M);
end

function M = calculate_mutual_inductance(spacing, length, width, wire_radius, num_turns, separation)
    % Initialize mutual inductance
    M = 0;
    
    % Calculate mutual inductance between each turn pair
    for i = 1:num_turns
        for j = 1:num_turns
            % Calculate dimensions of turn i in first coil
            l1 = length - 2*(i-1)*spacing;
            w1 = width - 2*(i-1)*spacing;
            
            % Calculate dimensions of turn j in second coil
            l2 = length - 2*(j-1)*spacing;
            w2 = width - 2*(j-1)*spacing;
            
            % Calculate mutual inductance between these turns using Neumann formula
            M = M + turn_mutual_inductance(l1, w1, l2, w2, separation);
        end
    end
end

function M = turn_mutual_inductance(l1, w1, l2, w2, separation)
    % Simplified mutual inductance calculation using Neumann formula
    % This is a basic approximation - for more accurate results, 
    % numerical integration should be used
    
    mu0 = 4*pi*1e-7; % Permeability of free space
    
    % Calculate geometric mean distance
    GMD = sqrt(separation^2 + (l1-l2)^2/12 + (w1-w2)^2/12);
    
    % Calculate mutual inductance
    M = mu0/(2*pi) * sqrt(l1*w1*l2*w2) / GMD;
end

function plot_coils(spacing, length, width, num_turns)
    figure;
    hold on;
    
    % Plot first coil
    z1 = 0;
    plot_single_coil(spacing, length, width, num_turns, z1, 'b');
    
    % Plot second coil
    z2 = 0.02; % 20mm separation
    plot_single_coil(spacing, length, width, num_turns, z2, 'r');
    
    view(3);
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Optimized Rectangular Coils');
    legend('Coil 1', 'Coil 2');
end

function plot_single_coil(spacing, length, width, num_turns, z, color)
    for i = 1:num_turns
        l = length - 2*(i-1)*spacing;
        w = width - 2*(i-1)*spacing;
        x = [-l/2, l/2, l/2, -l/2, -l/2];
        y = [-w/2, -w/2, w/2, w/2, -w/2];
        z_array = z * ones(size(x));
        plot3(x, y, z_array, color, 'LineWidth', 2);
    end
end 