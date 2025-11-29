function [result_table, losses_compar] = gwo_live_simulation()
% **Initialize Iteration Variables**
Pload_initial = evalin('base', 'no_drload_amplitude_24_1');  % [24x1] Initial load
convergence_threshold = 10;
iteration = 0;
max_iterations = 50;
% **GWO Parameters**
nVars = length(Pload_initial);
numWolves = 30;  
maxIter = 500;

%History
loss_history = [];
load_history = [];
best_loss = Inf;
best_load = [];

% **Define limits**
hourly_power_limit = [6;6;6;6;6;6;4;4;4;5;5;5;5;5;5;5;5;4;4;4;4;6;6;6] * 1000;
min_load = [2.0;2.0;2.0;2.0;2.0;2.0;5.5;5.5;5.5;4.5;4.5;4.5;4.5;4.5;4.5;4.5;4.5;5.5;5.5;5.5;5.5;2.0;2.0;2.0] * 1000;

if iteration == 0
% **Retrieve Initial Load from Workspace**

Pload_opt_live = Pload_initial;  % First time use the initial load
assignin('base', 'Pload_opt_live', Pload_opt_live); % Store it for Simulink

% Start the initial simulation
set_param('with_demand_res_final_v5_live', 'SaveOutput', 'on', 'ReturnWorkspaceOutputs', 'on');
fprintf('\n🔄 Running Initial Simulation...\n');
simOut = sim('with_demand_res_final_v5_live', 'ReturnWorkspaceOutputs', 'on');
assignin('base', 'with_demand_res_live', simOut);


% Ensure MATLAB waits for Simulink to finish
while ~evalin('base', 'exist("with_demand_res_live", "var")')
    pause(1);
end

    % Results of initial simulation (already have run)
    Pgrid_simulated = filter_values('Pgrid');
    Ploss_simulated = filter_values('Ploss');
    Pload_simulated = filter_values('Pload');

    % Initial losses
    current_loss = sum(Ploss_simulated);

    % ✅ **Store Results**
    loss_history = [loss_history; current_loss];
    load_history = [load_history; Pload_opt_live'];


end

original_total_energy = sum(Pload_initial);
iteration=1;

% Patience-based convergence
patience_counter = 0;
patience_limit = 5;  % Number of iterations allowed without improvement

% **Start Iterative Optimization**
while iteration < max_iterations

    iteration = iteration + 1;
    fprintf('📊 Iteration: %d\n', iteration);

    % Load optimization
    % Delete previous Pload_opt_live 
    if evalin('base', 'exist(''Pload_opt_live'', ''var'')')
    evalin('base', 'clear Pload_opt_live');
    end
    
    
    lb = max(Pload_initial - hourly_power_limit, min_load);
    ub = min(Pload_initial + hourly_power_limit, 20000);
    penalty_factor = 1 - exp(-sum(Pgrid_simulated) / sum(Pload_simulated));
    objectiveFunction = @(P_load) sum(Ploss_simulated) + ...
    (4000 * penalty_factor) * sqrt(abs(sum(Pload_simulated) - original_total_energy)) + ...
    (2500 * penalty_factor) * sum(max(0, Pload_simulated(2:6) - Pload_simulated(1:5)));
    
    % **Run Grey Wolf Optimization**
    Pload_opt_live = grey_wolf_optimizer(@(P) constraint_penalized_function(P, objectiveFunction), nVars, lb, ub, numWolves, maxIter);
    
   
    % **Ensure Pload_opt_live stays within bounds**
    Pload_opt_live = max(Pload_opt_live, lb);
    Pload_opt_live = min(Pload_opt_live, ub);
    Pload_opt_live = round(Pload_opt_live, -1);

    %  **Store Pload_opt_live for Simulink**
    assignin('base', 'Pload_opt_live', Pload_opt_live);
    
    % Delete previous results 
    if evalin('base', 'exist(''with_demand_res_live'', ''var'')')
    evalin('base', 'clear with_demand_res_live');
    end
   
    %  **Run new simulation Model with New Load
    set_param('with_demand_res_final_v5_live', 'SaveOutput', 'on', 'ReturnWorkspaceOutputs', 'on');
    fprintf('\n🔄 Running Simulation (Iteration %d)...\n', iteration);
    simOut = sim('with_demand_res_final_v5_live', 'ReturnWorkspaceOutputs', 'on');
    assignin('base', 'with_demand_res_live', simOut);
    %  Ensure MATLAB waits for new simulation to finish
    while ~evalin('base', 'exist("with_demand_res_live", "var")')
      pause(1);
    end

    %  New results in order to update the loss and load history
    Pgrid_simulated = filter_values('Pgrid');
    Ploss_simulated = filter_values('Ploss');
    Pload_simulated = filter_values('Pload');
    current_loss = sum(Ploss_simulated);

    % Update history
    loss_history = [loss_history; current_loss];
    load_history = [load_history; Pload_opt_live'];

    assignin('base', 'loss_history', loss_history);
    assignin('base', 'load_history', load_history);

     % ✅ Plot updated loss history
    figure(1);
    plot(1:length(loss_history), loss_history, 'b-o', 'LineWidth', 1.5);
    title('Live Optimization with Simulink');
    xlabel('Iteration');
    ylabel('Loss (W)');
    grid on;
    saveas(gcf, 'loss_plot.png');
    drawnow;

    % Update best solution & patience
    if current_loss < best_loss
       if current_loss > best_loss - convergence_threshold
           patience_counter = patience_counter + 1;
       else
           patience_counter = 0;
       end
       best_loss = current_loss;
       best_load = Pload_simulated;
    else
       patience_counter = 0;
    end

    
 
   % Check stopping condition
   if iteration > 10 && patience_counter >= patience_limit
        fprintf('✅ Patience window & convergence threshold met. Stopping optimization.\n');
        break;
   end

   if iteration == max_iterations
       fprintf('Patience window & convergence threshold DONT met. Stopping optimization.\n');
   end

end

% **Prepare Result Table**
result_table = table((1:nVars)', Pload_initial, best_load, ...
    'VariableNames', {'Step', 'Load_Before_W', 'Load_After_W'});

% **Compute Loss Comparison**
losses_compar = table(loss_history(1), best_loss, 'VariableNames', {'Initial_Losses', 'Final_Losses'});

% **Display Results**
disp(result_table);
disp(losses_compar);
fprintf('✅ Optimization process complete!\n');
end

%% **Constraint Penalty Function**
function penalty = constraint_penalized_function(P_load, objectiveFunction)
transition_penalty = sum(exp(max(0, abs(diff(P_load)) - 3.2)).^2) * 780;
penalty = objectiveFunction(P_load) + transition_penalty;
end

%% **Filtering Function**
function filtered_value = filter_values(signal_name)
    % Access the Simulink.SimulationOutput object
    simOut = evalin('base', 'with_demand_res_live');

    % Extract time and signal
    time = simOut.get('tout');
    signal = simOut.get(signal_name);

    % Reshape to column vectors
    time = reshape(time, [], 1);
    signal = reshape(signal, [], 1);

    % Initialize result
    num_points = 24;
    filtered_value = zeros(num_points, 1);

    % For each second i = 1 to 24, take values from (i-1)+0.85 to (i-1)+0.99
    for i = 1:num_points
        start_time = (i - 1) + 0.85;
        end_time = (i - 1) + 0.99;

        in_window = (time >= start_time) & (time < end_time);
        values = signal(in_window);

        if ~isempty(values)
            filtered_value(i) = mean(values, 'omitnan');
        else
            filtered_value(i) = mean(signal, 'omitnan');  % fallback if no data in window
        end
    end
end

%% **Grey Wolf Optimizer Function**
function [bestSol] = grey_wolf_optimizer(fitnessFunction, nVars, lb, ub, numWolves, maxIter)
wolves = repmat(lb', numWolves, 1) + rand(numWolves, nVars) .* repmat((ub - lb)', numWolves, 1);
scores = arrayfun(@(i) fitnessFunction(wolves(i, :)), 1:numWolves);
[~, idx] = sort(scores);
alpha = wolves(idx(1), :);

for iter = 1:maxIter
    a = 2 - iter * (2 / maxIter);
    for i = 1:numWolves
        for j = 1:nVars
            A1 = 2 * a * rand - a;
            C1 = 2 * rand;
            X1 = alpha(j) - A1 * abs(C1 * alpha(j) - wolves(i, j));
            wolves(i, j) = X1;
        end
    end
end
bestSol = alpha';
end