function [result_table] = gwo_nn_comb()
% Combined GWO and NN-based optimization for load adjustment based on available renewable energy

% Set random seed for consistent results
rng(1);

% Load filtered _no values from the workspace (convert from W to kW) & Ensure Column Vectors
Pwind_filtered_no = evalin('base', 'Pwind_filtered_no') / 1000; 
Ppv_filtered_no   = evalin('base', 'Ppv_filtered_no') / 1000;
Pgrid_filtered_no = evalin('base', 'Pgrid_filtered_no') / 1000; % Original Grid Power (No DR)
Pload_filtered_no = evalin('base', 'Pload_filtered_no') / 1000; % Original Load (No DR)

% Ensure all variables are Column Vectors (24×1)
Pwind_filtered_no = Pwind_filtered_no(:);
Ppv_filtered_no   = Ppv_filtered_no(:);
Pgrid_filtered_no = Pgrid_filtered_no(:);
Pload_filtered_no = Pload_filtered_no(:);

% Renewable Energy Calculation
P_renewable_no = Pwind_filtered_no + Ppv_filtered_no;

% Number of hours
num_hours = length(Pload_filtered_no);
original_total_energy = sum(Pload_filtered_no); % Keep total load demand stable

% GWO Parameters
nVars = num_hours;  
numWolves = 30;      % Number of search agents (wolves)
maxIter = 500;       % Maximum iterations

% Hourly flexibility limits (Ensure Column Vectors)
hourly_power_limit = [6; 6; 6; 6; 6; 6; ...
                      4; 4; 4; ...
                      5; 5; 5; 5; 5; 5; 5; 5; ...
                      4; 4; 4; 4; ...
                      6; 6; 6];

% Dynamic minimum load constraint (Ensure Column Vectors)
min_load = [2.0; 2.0; 2.0; 2.0; 2.0; 2.0; ...
            5.5; 5.5; 5.5; ...
            4.5; 4.5; 4.5; 4.5; 4.5; 4.5; 4.5; 4.5; ...
            5.5; 5.5; 5.5; 5.5; ...
            2.0; 2.0; 2.0];

% Bounds for load adjustment (Ensure Column Vectors)
lb = max(Pload_filtered_no - hourly_power_limit, min_load);
ub = Pload_filtered_no + hourly_power_limit;

% Ensure bounds correctness
if any(isnan(lb)) || any(isnan(ub))
    error('Error: lb or ub contains NaN values. Check constraints!');
end

% Adaptive Penalty Factor Based on Pre-Optimization Grid Power
penalty_factor = 1 - exp(-sum(Pgrid_filtered_no) / sum(Pload_filtered_no));

% Objective Function
objectiveFunction = @(P_load) sum(calculate_kgCO2(P_load, P_renewable_no)) + ...
   (4000 * penalty_factor) * sqrt(abs(sum(P_load) - original_total_energy)) + ...  
   (2500 * penalty_factor) * sum(max(0, P_load(2:6) - P_load(1:5)));

% Run GWO Optimization
[P_load_opt] = grey_wolf_optimizer(@(P) constraint_penalized_function(P, objectiveFunction), ...
                      nVars, lb, ub, numWolves, maxIter);

% Ensure `P_load_opt` is a column vector
P_load_opt = P_load_opt(:);

% Final bounds and rounding
P_load_opt = min(max(P_load_opt, lb), ub);
P_load_opt = round(P_load_opt, 3);

% Final energy correction
optimized_total_energy = sum(P_load_opt);
correction_factor = (original_total_energy / optimized_total_energy);
P_load_opt = round(P_load_opt * correction_factor, 3);

% Recalculate `P_grid_opt_no` AFTER GWO Optimization
P_grid_opt_no = P_load_opt - P_renewable_no; % Optimized Grid Power

% Store optimized load in workspace for Neural Network training
Pload_opt_gwo = P_load_opt;  % Assign the variable properly
assignin('base', 'Pload_opt_gwo', Pload_opt_gwo);

% Now proceed to NN optimization using the GWO results
% Neural Network-Based Load Optimization (NN-GWO)
Pwind = evalin('base', 'Pwind_filtered_no') / 1000; 
Ppv = evalin('base', 'Ppv_filtered_no') / 1000;
Pgrid_no = evalin('base', 'Pgrid_filtered_no') / 1000;
Pload_no = evalin('base', 'Pload_filtered_no') / 1000; 
Pload_opt_gwo = evalin('base', 'Pload_opt_gwo') / 1000; % From GWO Results

% Ensure Column Vectors
Pwind = Pwind(:);
Ppv = Ppv(:);
Pgrid_no = Pgrid_no(:);
Pload_no = Pload_no(:);
Pload_opt_gwo = Pload_opt_gwo(:);

% Feature Matrix (Inputs for NN)
X = [Pwind, Ppv, Pgrid_no, Pload_no]';

% Target Output (GWO Optimized Load)
Y = Pload_opt_gwo';

% Define Neural Network
numHiddenNeurons = 50; % Adjustable Hyperparameter
net = feedforwardnet(numHiddenNeurons);

% Set training parameters
net.trainParam.epochs = 500; % Set number of epochs
net.trainParam.max_fail = 20; % Allow more validation failures
net.trainParam.lr = 0.01; % Set learning rate
net.trainParam.showWindow = true; % Show training window
net.trainParam.goal = 0.001; % Goal for the training error
net.trainParam.min_grad = 1e-6; % Minimum gradient to stop training

% Train Neural Network
net.trainFcn = 'trainscg'; % Using Scaled Conjugate Gradient
net = train(net, X, Y);

% Predict Optimized Load Using NN
Pload_opt_nn = net(X)';

% Ensure Column Vector
Pload_opt_nn = Pload_opt_nn(:);

% Improved Fix for Exact Energy Balance
original_total_energy = sum(Pload_no);
optimized_total_energy = sum(Pload_opt_nn);

while abs(optimized_total_energy - original_total_energy) > 0.001
    correction_factor = original_total_energy / optimized_total_energy;
    Pload_opt_nn = round(Pload_opt_nn * correction_factor, 3);
    optimized_total_energy = sum(Pload_opt_nn);
end

% Ensure Predicted Load Stays Within Bounds
hourly_power_limit = [6; 6; 6; 6; 6; 6; ...
                      4; 4; 4; ...
                      5; 5; 5; 5; 5; 5; 5; 5; ...
                      4; 4; 4; 4; ...
                      6; 6; 6];

min_load = [2.0; 2.0; 2.0; 2.0; 2.0; 2.0; ...
            5.5; 5.5; 5.5; ...
            4.5; 4.5; 4.5; 4.5; 4.5; 4.5; 4.5; 4.5; ...
            5.5; 5.5; 5.5; 5.5; ...
            2.0; 2.0; 2.0];

lb = max(Pload_no - hourly_power_limit, min_load);
ub = Pload_no + hourly_power_limit;
Pload_opt_nn = min(max(Pload_opt_nn, lb), ub);

% Recalculate Grid Power After NN Optimization
Pgrid_opt_nn = Pload_opt_nn - (Pwind + Ppv);

% Store NN-optimized load in the workspace
Pload_opt_nngwo_comb_co2 = Pload_opt_nn;  % Assign the variable properly
assignin('base', 'Pload_opt_nngwo_comb_co2', Pload_opt_nngwo_comb_co2);

% Cost & CO2 Calculations
cost_before = calculate_hourly_cost(Pload_no, Pgrid_no);
cost_after = calculate_hourly_cost(Pload_opt_nn, Pgrid_opt_nn);

kgCO2_before = calculate_kgCO2(Pload_no, Pgrid_no);
kgCO2_after = calculate_kgCO2(Pload_opt_nn, Pgrid_opt_nn);

% Final Energy Summary
total_row = [25, sum(Pload_no), sum(Pload_opt_nn), sum(cost_before), sum(cost_after), sum(kgCO2_before), sum(kgCO2_after)];

% Create Result Table
result_table = array2table([[(1:length(Pload_no))', Pload_no, Pload_opt_nn, cost_before, cost_after, kgCO2_before, kgCO2_after]; total_row], ...
    'VariableNames', {'Hour', 'Load_Before_kW', 'Load_After_NN_kW', 'Cost_Before_USD', 'Cost_After_USD', 'kgCO2_Before', 'kgCO2_After'});

% Display Results
disp(result_table);
fprintf('\nOriginal Total Load Energy: %.3f kWh\n', sum(Pload_no));
fprintf('NN Optimized Total Load Energy: %.3f kWh\n', sum(Pload_opt_nn));

end

%% Grey Wolf Optimizer Function
function [bestSol] = grey_wolf_optimizer(fitnessFunction, nVars, lb, ub, numWolves, maxIter)
    % Initialize wolf pack positions within bounds
    wolves = repmat(lb', numWolves, 1) + rand(numWolves, nVars) .* repmat((ub - lb)', numWolves, 1);

    % Evaluate initial wolves
    scores = arrayfun(@(i) fitnessFunction(wolves(i, :)), 1:numWolves);

    % Identify alpha, beta, delta wolves
    [~, idx] = sort(scores);
    alpha = wolves(idx(1), :);
    beta = wolves(idx(2), :);
    delta = wolves(idx(3), :);

    for iter = 1:maxIter
        a = 2 - iter * (2 / maxIter); % Decreasing exploration factor

        for i = 1:numWolves
            for j = 1:nVars
                A1 = 2 * a * rand - a;
                C1 = 2 * rand;
                X1 = alpha(j) - A1 * abs(C1 * alpha(j) - wolves(i, j));

                A2 = 2 * a * rand - a;
                C2 = 2 * rand;
                X2 = beta(j) - A2 * abs(C2 * beta(j) - wolves(i, j));

                A3 = 2 * a * rand - a;
                C3 = 2 * rand;
                X3 = delta(j) - A3 * abs(C3 * delta(j) - wolves(i, j));

                wolves(i, j) = (X1 + X2 + X3) / 3;
            end
        end
    end

    bestSol = alpha';
end

%% Cost Calculation Functions
function cost = calculate_hourly_cost(load, P_grid)
    t = sum(P_grid > 0);
    P_grid_max = max(P_grid);
    cost = arrayfun(@(P) compute_cost(P, t, P_grid_max), P_grid);
end
% **Constraint Penalty Function for GWO**
function penalty = constraint_penalized_function(P_load, objectiveFunction)
    % This penalty function adds a penalty for violating constraints such as large transitions in load.
    
    % Transition penalty for load variations
    transition_penalty = sum(exp(max(0, abs(diff(P_load)) - 3.2)).^2) * 780;
    
    % Calculate the objective function value (from cost and other factors)
    penalty = objectiveFunction(P_load) + transition_penalty;
end


function cost = compute_cost(P_grid, t, P_grid_max)
    Cinv = 1.78;
    Cmnt = ((1.7 * (12 * t / 4)) + (24.1 * P_grid_max)) / (365 * 24);
    if P_grid > 0
        Cop = (0.51 * P_grid + 1.87) + (0.04 * P_grid + 0.15) + (4.1e-5 * P_grid) + (0.1e-5 * P_grid);
        cost = Cinv + Cmnt + Cop;
    else
        cost = Cinv + Cmnt;
    end
end

%% CO2 Calculation Function
function kgCO2 = calculate_kgCO2(load, P_grid)
    kgCO2 = arrayfun(@(P) (P > 0) * (2.51 * (0.26 * P + 0.95)), P_grid);
end
