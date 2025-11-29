function [Pwind_filtered_with_loss, Ppv_filtered_with_loss, Pgrid_filtered_with_loss, Pload_filtered_with_loss, Psupply_filtered_with_loss, Pin_trans_filtered_with_loss, Pout_trans_filtered_with_loss, Ploss_filtered_with_loss, Ploss_LVline_filtered_with_loss, Ploss_MVline_filtered_with_loss, Ploss_trans_filtered_with_loss, ...
          Qwind_filtered_with_loss, Qpv_filtered_with_loss, Qgrid_filtered_with_loss, Qload_filtered_with_loss, Qsupply_filtered_with_loss, Qin_trans_filtered_with_loss, Qout_trans_filtered_with_loss, Qloss_filtered_with_loss, Qloss_LVline_filtered_with_loss, Qloss_MVline_filtered_with_loss, Qloss_trans_filtered_with_loss,...
          Iwind_filtered_with_loss, Ipv_filtered_with_loss, Igrid_filtered_with_loss, Iload_filtered_with_loss, Isupply_filtered_with_loss, Iin_trans_rms_filtered_with_loss, Iout_trans_rms_filtered_with_loss, ...
          f_wind_filtered_with_loss, f_pv_filtered_with_loss, f_grid_filtered_with_loss, f_load_filtered_with_loss, f_supply_filtered_with_loss, ...
          Vwind_filtered_with_loss, Vpv_filtered_with_loss, Vgrid_filtered_with_loss, Vload_filtered_with_loss, V_supply_filtered_with_loss, Vin_trans_filtered_with_loss, Vout_trans_filtered_with_loss, ...
          time_filterd_with_loss] ...
           = filter_signals(with_demand_res_loss)

    % If 'with_demand_res_loss' is witht provided, retrieve it from the workspace
    if nargin < 1
        if evalin('base', 'exist(''with_demand_res_loss'', ''var'')')
            with_demand_res_loss = evalin('base', 'with_demand_res_loss');
        else
            error('The variable "with_demand_res_loss" is missing. Run the Simulink simulation first.');
        end
    end

    % Extract variables from Simulink output structure
    variable_names = {'Pwind', 'Ppv', 'Pgrid', 'Pload', 'Psupply', 'Pin_trans', 'Pout_trans', 'Ploss', 'Ploss_LVline', 'Ploss_MVline', 'Ploss_trans', ...
                      'Qwind', 'Qpv', 'Qgrid', 'Qload', 'Qsupply', 'Qin_trans', 'Qout_trans', 'Qloss', 'Qloss_LVline', 'Qloss_MVline', 'Qloss_trans', ...
                      'Iwind', 'Ipv', 'Igrid', 'Iload', 'Isupply', 'Iin_trans_rms', 'Iout_trans_rms' ...
                      'f_wind', 'f_pv', 'f_grid', 'f_load', 'f_supply', ...
                      'Vwind', 'Vpv', 'Vgrid', 'Vload', 'V_supply', 'Vin_trans', 'Vout_trans'};
                  
    % Create "_with" versions of all variables
    for i = 1:length(variable_names)
        var_name_with_loss = [variable_names{i}, '_with_loss'];
        eval([var_name_with_loss ' = with_demand_res_loss.get(''' variable_names{i} ''');']);
    end

    time_with_loss = with_demand_res_loss.get('time');

    % Ensure all variables are column vectors (24×1)
    for i = 1:length(variable_names)
        var_name_with_loss = [variable_names{i}, '_with_loss'];
        eval([var_name_with_loss ' = reshape(' var_name_with_loss ', [], 1);']);
    end

    time_with_loss = reshape(time_with_loss, [], 1);

    % Initialize filtered outputs as column vectors (24×1)
    num_hours = 24;
    time_filtered_with_loss = zeros(num_hours, 1);
    for i = 1:length(variable_names)
        var_name_filtered_with_loss = [variable_names{i}, '_filtered_with_loss'];
        eval([var_name_filtered_with_loss ' = zeros(num_hours, 1);']);
    end

    % Process each hour separately (1–24)
    for i = 1:num_hours  
        % Find indices where time step is within the correct range
        valid_indices = (mod(time_with_loss, 1) >= 0.85) & (mod(time_with_loss, 1) < 0.99) & (time_with_loss >= i-1) & (time_with_loss < i);

        % Compute averages for all extracted variables
        for j = 1:length(variable_names)
            var_name_with_loss = [variable_names{j}, '_with_loss'];
            var_name_filtered_with_loss = [variable_names{j}, '_filtered_with_loss'];
            eval([var_name_filtered_with_loss '(i, 1) = mean(' var_name_with_loss '(valid_indices), ''omitnan'');']);
        end
        if any(valid_indices)
        time_filtered_with_loss(i, 1) =round(time_with_loss(find(valid_indices, 1, 'last')));
    else
        time_filtered_with_loss(i, 1) = NaN;  % Handle cases where with valid time is found
    end
    end

    % Assign both raw and filtered variables to the workspace as **column vectors**
    for i = 1:length(variable_names)
        var_name_with_loss = [variable_names{i}, '_with_loss'];
        var_name_filtered_with_loss = [variable_names{i}, '_filtered_with_loss'];
        assignin('base', var_name_with_loss, eval(var_name_with_loss));  % Raw data
        assignin('base', var_name_filtered_with_loss, eval(var_name_filtered_with_loss));  % Filtered data
        assignin('base', 'time_filtered_with_loss', time_filtered_with_loss);
    end

    assignin('base', 'time_with_loss', time_with_loss);

    disp('Filtered values with "_with_loss" and "_filtered_with_loss" suffixes have been saved to the workspace as column vectors.');
    disp('Filtered time (time_filtered_with_loss) saved to the workspace.');
end
