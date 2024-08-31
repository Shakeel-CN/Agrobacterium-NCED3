function bac_behavior_nced_8_14
    % Model parameters
    params.L = 10;           % Space dimensions (L x L)
    params.dx = 0.1;         % Spatial step size
    params.dy = params.dx;   % Spatial step size (same as dx)
    params.dt = 0.01;        % Time step
    params.t_end = 100;      % Total simulation time

    % Diffusion coefficients
    params.D_nced = 0.1;     % Diffusion coefficient of NCED
    params.D_n = 0.05;       % Diffusion coefficient of nutrient

    % Rate constants
    params.k_p = 0.5;        % Rate constant for NCED production
    params.gamma_nced = 0.05;% Rate constant for NCED degradation
    params.k_n = 0.2;        % Rate constant for nutrient consumption
    params.rho = 0.1;        % Rate constant for bacterial death
    params.k_aba = 0.1;      % Rate constant for ABA synthesis from NCED

    % Bacterial movement parameters
    params.v = 10e-6;        % Average speed of bacteria
    params.sigma = 0.5;      % Standard deviation of random walk step size

    % Define chemotaxis function (example)
    chi = @(x,y) exp(-(x-params.L/2).^2 / (0.5*params.L)^2) .* exp(-(y-params.L/2).^2 / (0.5*params.L)^2); % Creates a peak at center

    % Initialize spatial grids
    [X,Y] = meshgrid(0:params.dx:params.L, 0:params.dy:params.L);

    % Initialize variables
    NCED = 0.1*ones(size(X)); % NCED concentration
    N = 2*ones(size(X));      % Nutrient concentration
    B = 0.1*ones(size(X));    % Bacterial density
    ABA = zeros(size(X));     % ABA concentration (initially zero)

    % Initial bacterial positions (random)
    B_x = rand(1000,1) * params.L;
    B_y = rand(1000,1) * params.L;

    % Pre-calculate Laplacian operator
    laplacian_operator = del2(eye(size(X)), params.dx) + del2(eye(size(Y)), params.dy);

   % Simulation loop
    for t = 0:params.dt:params.t_end
        % Solve diffusion equations for NCED and nutrient
        NCED = NCED + params.dt * (params.D_nced * laplacian_operator * NCED + params.k_p * B .* min(N, 1) - params.gamma_nced * NCED);
        N = N + params.dt * (params.D_n * laplacian_operator * N - params.k_n * B .* N .* (N > 0.1));

        % Update bacterial positions (biased movement)
        dx = randn(1000,1) * params.sigma;
        dy = randn(1000,1) * params.sigma;
        dx = dx + params.v * params.dt * chi(B_x, B_y) .* cos(atan2(dy, dx));
        dy = dy + params.v * params.dt * chi(B_x, B_y) .* sin(atan2(dy, dx));
        B_x = B_x + dx;
        B_y = B_y + dy;
        B_x = mod(B_x, params.L); % Wrap around edges
        B_y = mod(B_y, params.L);

        % Update bacterial density
        B = update_bacterial_density(B, B_x, B_y, params, laplacian_operator, N);

        % Update ABA concentration based on NCED presence
        ABA = ABA + params.dt * (params.k_aba * NCED - params.gamma_nced * ABA);

        % Visualize results in separate figures
        visualize_results(X, Y, NCED, B, N, B_x, B_y, ABA, params, t, 'data.csv');

        pause(0.1); % Optional pause for slower visualization
    end
end

% Define function to update bacterial density
function B = update_bacterial_density(B, B_x, B_y, params, laplacian_operator, N)
    for i = 1:numel(B_x)
        % Calculate grid indices (corrected clamping)
        x_idx = max(1, min(floor(B_x(i)/params.dx) + 1, size(B,1)));
        y_idx = max(1, min(floor(B_y(i)/params.dy) + 1, size(B,2)));

        % Increase density at the bacteria's position
        B(x_idx, y_idx) = B(x_idx, y_idx) + 1;
    end
    
    % Update bacterial density based on movement, production, and death
    B = B + params.dt * (params.D_n * laplacian_operator * B - params.k_n * B .^2 .* (N > 0.1) - params.rho * B);
end

function visualize_results(X, Y, NCED, B, N, B_x, B_y, ABA, params, t, filename)
    figure;
    
    subplot(2,3,1);
    surf(X, Y, NCED);
     shading interp;
    title('NCED Concentration');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    zlabel('Concentration');
    colormap(parula); % lighter colormap
    
    subplot(2,3,2);
    surf(X, Y, B);
     shading interp;
    title('Bacterial Density');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    zlabel('Density');
    colormap(parula); % lighter colormap
    
    subplot(2,3,3);
    surf(X, Y, N);
     shading interp;
    title('Nutrient Availability');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    zlabel('Availability');
    colormap(parula); % lighter colormap
    
    subplot(2,3,4);
    surf(X, Y, ABA);
     shading interp;
    title('ABA Concentration');
    xlabel('X (nm)');
    ylabel('Y (nm)');
    zlabel('Concentration');
    colormap(parula); % lighter colormap
    
    % Plot bacterial trajectories
    hold on;
    scatter3(B_x(:), B_y(:), zeros(size(B_x(:))), 'r.'); % Plot initial positions
    for i = 1:size(B_x, 1)
        plot3(B_x(i, :), B_y(i, :), zeros(size(B_x(i, :))), 'b-'); % Plot trajectories
    end
    hold off;
    
    % Export data to CSV file
    data = [NCED(:), B(:), N(:), ABA(:)];
    writematrix(data,filename);
    
    pause(0.5);
end
