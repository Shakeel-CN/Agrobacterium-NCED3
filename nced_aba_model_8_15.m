function nced_aba_model_8_15()
    % Set parameters for agent-based model
    time_constant = 500; % seconds
    step_length = 5; % microns
    diffusion_coefficient = 10e-3; % um^2/s (can be calculated based on viscosity and temperature)
    initial_bacteria_concentration = 1e8; % CFU/mL
    absorption_rate = 0.005; % rate at which bacteria is absorbed (dimensionless)
    gain_Tx = 1.5;  % dimensionless
    gain_B = 1;  % dimensionless (ensure realistic impact)
    gain_Rx = 0.9;  % dimensionless
    delay_Tx = 0.1; % seconds
    delay_B = 0.005; % seconds
    delay_Rx = 0.9; % seconds
    receiver_location = 15; % microns
    
    % Parameters for the differential equations
    k_u = 0.05; % Gene uptake rate constant (1/s)
    k_d = 0.01; % Gene degradation rate constant (1/s)
    k_e = 0.1; % Enzyme production rate constant (1/s)
    k_m = 0.05; % Enzyme degradation rate constant (1/s)
    k_ABA = 0.02; % ABA production rate constant (1/s)
    k_ABA_d = 0.01; % ABA degradation rate constant (1/s)
    
    % Simulate agent diffusion
    time_steps = 1000;
    agent_position = simulateAgentDiffusion(time_steps, step_length);

    % Adjust time vector to match the number of time steps
    t_sim = linspace(0, time_constant, time_steps);

    % Simulate signal emission (gene delivery)
    signal_Tx = simulateGeneDelivery(agent_position, gain_Tx, t_sim);

    % Simulate signal propagation (bacterial movement and efficiency)
    signal_B = simulateGenePropagation(signal_Tx, diffusion_coefficient, time_steps, gain_B, initial_bacteria_concentration, absorption_rate);

    % Simulate signal reception (gene uptake by plant cells)
    [signal_Rx, gene_concentration, enzyme_concentration, aba_production] = simulateGeneUptake(signal_B, agent_position, receiver_location, gain_Rx, k_u, k_d, k_e, k_m, k_ABA, k_ABA_d, t_sim);

    % t_emit and t_receive should be determined dynamically based on the signal
    [t_emit, t_receive] = determineEmitReceiveTimes(signal_Tx, signal_Rx, t_sim);

    % Calculate propagation latency
    propagation_latency = t_receive - t_emit;

    % Calculate total latency
    total_latency = delay_Tx + propagation_latency + delay_Rx;

    fprintf('Total Latency: %.2f seconds\n', total_latency);

    % Plot relevant data with colorful lines
    figure;
    subplot(6, 1, 1);
    plot(t_sim, agent_position, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Bacterial Position (\mum)');
    title('(a) Agrobacterium Position');
    grid on;

    subplot(6, 1, 2);
    plot(t_sim, signal_Tx, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Gene Delivery Signal (ng/\muL)');
    title('(b) NCED3 Gene Delivery by Bacteria');
    hold on;
    plot([t_emit t_emit], [min(signal_Tx) max(signal_Tx)], 'k--');
    hold off;
    grid on;

    subplot(6, 1, 3);
    plot(t_sim, signal_B, 'g', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Gene Propagation Signal (cm^2/h)');
    title('(c) NCED3 Gene Propagation through Bacterial Movement');
    hold on;
    plot([t_receive t_receive], [min(signal_B) max(signal_B)], 'k--');
    hold off;
    grid on;

    subplot(6, 1, 4);
    plot(t_sim, gene_concentration, 'b', 'LineWidth', 2);
   xlabel('Time (s)');
    ylabel('Gene Concentration (ng/\muL)');
    title('(d) Concentration of NCED3 Gene Over Time');
    grid on;

    subplot(6, 1, 5);
    plot(t_sim, signal_Rx, 'm', 'LineWidth', 2);
   xlabel('Time (s)');
    ylabel('Gene Uptake Signal (copies/cell)');
    title('(e) NCED3 Gene Uptake by Plant Cells');
    hold on;
    plot([t_receive t_receive], [min(signal_Rx) max(signal_Rx)], 'k--');
    hold off;
    grid on;

    subplot(6, 1, 6);
    plot(t_sim, aba_production, 'c', 'LineWidth', 2);
     xlabel('Time (s)');
    ylabel('ABA Production (\mumol)');
    title('(f) ABA Production Over Time');
    grid on;

    % Agent Diffusion 
    function agent_position = simulateAgentDiffusion(time_steps, step_length)
        agent_position = zeros(1, time_steps);
        for i = 2:time_steps
            % Motion equation for agent
            agent_position(i) = agent_position(i-1) + step_length * (-1 + 2 * rand);
        end
    end

    function signal_Tx = simulateGeneDelivery(agent_position, gain_Tx, t_sim)
        % Simulate gene delivery based on the spatial concentration ratio
        signal_Tx = zeros(1, length(agent_position));
        max_signal = gain_Tx * max(abs(agent_position));
        for i = 1:length(agent_position)
            if agent_position(i) > 0
                % Create initial low, then high, then low pattern
                signal_Tx(i) = gain_Tx * abs(agent_position(i)) * sin(pi * t_sim(i) / max(t_sim));
            else
                % Create initial low, then high, then low pattern
                signal_Tx(i) = gain_Tx * abs(agent_position(i)) * 0.5 * sin(pi * t_sim(i) / max(t_sim));
            end
        end
        signal_Tx = max_signal * signal_Tx / max(signal_Tx); % Normalize to max signal
    end

    function signal_B = simulateGenePropagation(signal_Tx, diffusion_coefficient, time_steps, gain_B, initial_bacteria_concentration, absorption_rate)
        % Simulate gene propagation using a simplified model
        signal_B = zeros(1, length(signal_Tx));
        for i = 2:length(signal_Tx)
            % Propagation calculation
            signal_B(i) = gain_B * diffusion_coefficient * signal_Tx(i); % Convert to appropriate units
            % Apply absorption
            signal_B(i) = signal_B(i) * (1 - absorption_rate);
        end
    end

    function [signal_Rx, gene_concentration, enzyme_concentration, aba_production] = simulateGeneUptake(signal_B, agent_position, receiver_location, gain_Rx, k_u, k_d, k_e, k_m, k_ABA, k_ABA_d, t_sim)
        % Simulate gene uptake and ABA production
        signal_Rx = zeros(1, length(agent_position));
        gene_concentration = zeros(1, length(agent_position));
        enzyme_concentration = zeros(1, length(agent_position));
        aba_production = zeros(1, length(agent_position));

        for i = 2:length(agent_position)
            if agent_position(i) >= receiver_location
                signal_Rx(i) = gain_Rx * signal_B(i); % Maintain appropriate units
            else
                signal_Rx(i) = gain_Rx * signal_B(i) * 0.5; % Maintain appropriate units
            end

            % Gene uptake and degradation
            gene_concentration(i) = gene_concentration(i-1) + k_u * signal_B(i) - k_d * gene_concentration(i-1);

            % Enzyme production and degradation
            enzyme_concentration(i) = enzyme_concentration(i-1) + k_e * gene_concentration(i) - k_m * enzyme_concentration(i-1);

            % ABA production and degradation
            aba_production(i) = aba_production(i-1) + k_ABA * enzyme_concentration(i) - k_ABA_d * aba_production(i-1);
        end
    end

    function [t_emit, t_receive] = determineEmitReceiveTimes(signal_Tx, signal_Rx, t_sim)
        % Dynamic determination of t_emit and t_receive based on signal data
        t_emit = t_sim(find(signal_Tx > 0, 1));
        t_receive = t_sim(find(signal_Rx > 0, 1, 'last'));
    end
end
