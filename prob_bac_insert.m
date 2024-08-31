% Stomatal Density (number of stomata per unit leaf area)
stomatal_density = 300; % Example value, adjust as needed

% Agrobacterium Arrival Rate (number of Agrobacterium arriving per unit time)
lambda = 50; % can be adjust as needed

% Agrobacterium Speed (μm/s)
agrobacterium_speed = 10; % Can be adjust as needed

% Agrobacterium Diffusion Coefficient (μm^2/s)
agrobacterium_diffusion_coeff = 20; % Can adjust as needed

% Viscosity of Leaf Surface (Pa·s)
viscosity = 0.01; % Can be adjust as needed

% Stomatal Size (μm)
stomatal_size = 30; % Can be adjust as needed

% Stomatal Accessibility Probability
stomatal_accessibility_prob = 0.2; % Can be adjust as needed

% Time Interval
time_interval = 1; % Values (in hours), Can be adjust as needed

% Number of Monte Carlo Simulations
num_simulations = 1000; % Adjust as needed
% Define the leaf area (e.g., 10 mm x 10 mm)
leaf_area = [2, 2]; % Can be as needed

% Generate random stomatal positions within the leaf area
num_stomata = poissrnd(stomatal_density * prod(leaf_area));
stomatal_positions = rand(num_stomata, 2) .* repmat(leaf_area, num_stomata, 1);
% Initialize variables
num_internalized = 0;
droplet_positions = [];

% Perform Monte Carlo Simulations
for i = 1:num_simulations
    % Generate random droplet position
    droplet_position = rand(1, 2) .* leaf_area;
    droplet_positions = [droplet_positions; droplet_position];

    % Check if the droplet lands on a stoma
    distances = sqrt(sum(bsxfun(@minus, stomatal_positions, droplet_position).^2, 2));
    if any(distances < stomatal_size / 2) % Check if the droplet lands on a stoma
        % Calculate the Agrobacterium arrival rate
        agrobacterium_arrival = poissrnd(lambda / time_interval);

        % Simulate Agrobacterium movement and internalization
        for j = 1:agrobacterium_arrival
            % Calculate the Agrobacterium displacement due to diffusion and speed
            dx = sqrt(2 * agrobacterium_diffusion_coeff * time_interval / viscosity) * randn();
            dy = sqrt(2 * agrobacterium_diffusion_coeff * time_interval / viscosity) * randn();
            droplet_position = droplet_position + [agrobacterium_speed * time_interval, 0] + [dx, dy];

            % Check if the Agrobacterium has entered the stoma
            if any(sqrt(sum(bsxfun(@minus, stomatal_positions, droplet_position).^2, 2)) < stomatal_size / 2)
                num_internalized = num_internalized + 1;
                break; % Exit the inner loop once internalization occurs
            end
        end
    end
end
% Calculate the likelihood of successful Agrobacterium internalization
likelihood_of_internalization = num_internalized / num_simulations;

% Display the results
fprintf('Likelihood of successful Agrobacterium internalization: %.4f\n', likelihood_of_internalization);

% Optionally, you can visualize the droplet positions and internalized stomata
figure;
plot(stomatal_positions(:, 1), stomatal_positions(:, 2), 'ro', 'MarkerSize', 4);
hold on;
plot(droplet_positions(:, 1), droplet_positions(:, 2), 'b.');
title('Leaf Area with Stomatal and Droplet Positions');
xlabel('Length (cm)');
ylabel('Width (cm)');
legend('Stomata', 'Droplet Positions');
