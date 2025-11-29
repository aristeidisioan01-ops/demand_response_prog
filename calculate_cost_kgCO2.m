function [result_table, total_cost, total_co2] = calculate_cost_co2()
    % Retrieve Pgrid from workspace and ensure it's in kW

    Pgrid = evalin('base', 'Pgrid_filtered_with_loss') / 1000;
    Pgrid = Pgrid(:);  % Ensure column vector

    % Initialization
    hours = numel(Pgrid);
    t = sum(Pgrid > 0);         % Time steps where grid is used
    Pgrid_max = max(Pgrid);     % Max grid power for maintenance cost

    % Preallocate
    cost = zeros(hours, 1);
    kgCO2 = zeros(hours, 1);

    % Calculate cost & CO2 for each hour
    for i = 1:hours
        P = Pgrid(i);
        Cinv = 1.78;
        Cmnt = ((1.7 * (12 * t / 4)) + (24.1 * Pgrid_max)) / (365 * 24);

        if P > 0
            Cop = (0.51 * P + 1.87) + (0.04 * P + 0.15) + (4.1e-5 * P) + (0.1e-5 * P);
            cost(i) = Cinv + Cmnt + Cop;
            kgCO2(i) = 2.51 * (0.26 * P + 0.95);
        else
            cost(i) = Cinv + Cmnt;
            kgCO2(i) = 0;
        end
    end

    % Total values
    total_cost = sum(cost);
    total_co2 = sum(kgCO2);

    % Assign to workspace (optional)
    assignin('base', 'cost_hourly', cost);
    assignin('base', 'kgCO2_hourly', kgCO2);
    assignin('base', 'total_cost', total_cost);
    assignin('base', 'total_kgCO2', total_co2);

    % Table
    result_table = table((1:hours)', Pgrid, cost, kgCO2, ...
        'VariableNames', {'Hour', 'Pgrid_kW', 'Cost_USD', 'kgCO2'});

    % Display results
    disp(result_table);
    fprintf('\n💸 Total Cost: %.3f USD\n', total_cost);
    fprintf('🌍 Total CO₂ Emissions: %.3f kgCO₂\n', total_co2);
end
