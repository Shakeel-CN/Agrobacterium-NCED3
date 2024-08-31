% Define the domain and initial conditions
A = 0.001; % Area of domain
N = 100; % Number of grid points
dx = A/N; % Grid spacing
x = linspace(0, A, N+1); % Grid points

% Varying parameters
D = 50e-12; % Bacterial diffusion constant (average value)
velocity = 10e-6; % Velocity of E.coli in m/s
pore_size_range = linspace(10e-6, 80e-6, 100); % Pore size range in microm
pore_length_range = linspace(7.6e-6, 26.914e-6, 100); % Pore length range in microm
temp_range_C = linspace(20, 32, 100); % Temperature range in Celsius
temp_range_K = temp_range_C + 273.15; % Convert temperature to Kelvin
time_range = linspace(0, 100, 100); % Time range in seconds
% Preallocate matrices for time, growth rate, and diffusion time
pore_matrix = zeros(length(pore_size_range),length(temp_range_K));
diffusion_time_matrix = zeros(length(temp_range_K), length(pore_size_range));

% Preallocate matrices for total time and concentration
total_time_matrix = zeros(length(pore_size_range), length(pore_length_range), length(temp_range_K));
concentration_matrix = zeros(N+1, length(pore_size_range), length(pore_length_range), length(temp_range_K));

% Free water diffusion coefficient for agrobacterium at 25 degrees Celsius
D_free = 50e-6; % m^2/s

% Loop through all parameter combinations
for p_idx = 1:length(pore_size_range)
    for l_idx = 1:length(pore_length_range)
        for temp_idx = 1:length(temp_range_K)
            pore_size = pore_size_range(p_idx);
            pore_length = pore_length_range(l_idx);
            T = temp_range_K(temp_idx);
            temp_effect = exp(-(T - min(temp_range_K))/(max(temp_range_K) - min(temp_range_K))); % Temperature effect
            
            % Calculate the specific growth rate of Agrobacterium
            growth_rate = 0.27/3600; % Set Avg growth rate to 0.75 h^-1 and change from hr to sec

            % Calculate the adjusted velocity based on pore size and temperature
            adjusted_velocity = velocity * (pore_size / max(pore_size_range)) * temp_effect;

            % Compute the effective diffusion coefficient using given equation
            D_eff = D_free;

            % Calculate time taken for bacteria to pass through porous medium
            time_diffusion = (A^2) / (2 * D_eff);
            time_advection = A / adjusted_velocity;
            
            % Include running and tumbling times
            running_time = 1.25; % sec
            tumbling_time = 0.17; % sec
            
            time_total = time_diffusion + time_advection + running_time - tumbling_time;

            % Calculate the concentration profile using the diffusion equation
            concentration = zeros(N+1, 1);
            for i = 1:N+1
                concentration(i) = exp(-((x(i) - A/2)^2) / (4 * D_eff * time_total)) / sqrt(4 * pi * D_eff * time_total);
            end

            % Update the concentration profile based on growth
            growth_factor = exp(growth_rate * time_total);
            concentration = concentration * growth_factor;

            % Store the results in the matrices
            total_time_matrix(p_idx, l_idx, temp_idx) = time_total;
            concentration_matrix(:, p_idx, l_idx, temp_idx) = concentration;
        end
    end
end

% Calculate the average total time
average_total_time = mean(total_time_matrix(:));

% Create meshgrids for visualization
[Pore_grid, Temp_grid] = meshgrid(pore_size_range, temp_range_C);

% Fixing pore length for plotting
pore_length_idx = 1;

% Extracting the data for fixed pore length
average_total_time_fixed_length = squeeze(total_time_matrix(:, pore_length_idx, :));
% Display the average total time
fprintf('The average total propaagation time taken by Agrobacterium to pass through the porous medium is %.2f seconds.\n', average_total_time);

% Plotting the 3D surface
figure;
surf(temp_range_C, pore_size_range , average_total_time_fixed_length);
xlabel('Temperature (C)');
ylabel('Pore Size (\mum)');
zlabel('Average Total Time (s)');
title('Effect of Temperature on Average Total Propagation Time for Agrobacterium');

% Create a 3D surface plot for time taken by Agrobacterium
[X, Y] = meshgrid(pore_size_range, pore_length_range);
figure('Color', 'white'); % Setting the background to white
h = surf(X, Y, squeeze(total_time_matrix(:, :, 1)), 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Improve lighting
light('Position', [-1 0 1], 'Style', 'infinite');
shading interp;
material dull;

% Add labels and title
xlabel('Pore Size (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pore Length (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Average Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Average Total Propagation Time for Agrobacterium to Pass Through Stomata of Wheat Leaf', 'FontSize', 14, 'FontWeight', 'bold');

% Set colormap
colormap(parula(256));

% Adjust view angle for better perspective
view(3);
% Create a 3D surface plot for time taken by Agrobacterium
[X, Y] = meshgrid(pore_size_range, pore_length_range);
figure;
surf(X, Y, squeeze(total_time_matrix(:, :, 1))); % Choose a specific temperature slice

% Add labels and title
xlabel('Pore Size (\mum)');
ylabel('Pore Length (\mum)');
zlabel('Average Total Time (s)');
title('Average Total Propagation Time for Agrobacterium to Pass Through Stomata of Wheat Leaf', 'FontSize', 14, 'FontWeight', 'bold');

% Customize axis and grid
axis tight;
grid on;
set(gca, 'GridAlpha', 0.3, 'Box', 'off', 'FontSize', 12, 'LineWidth', 1.2);

