function [Pwind_filtered_with, Ppv_filtered_with, Pgrid_filtered_with, Pload_filtered_with, Psupply_filtered_with, Pin_trans_filtered_with, Pout_trans_filtered_with, Ploss_filtered_with, Ploss_LVline_filtered_with, Ploss_MVline_filtered_with, Ploss_trans_filtered_with, ...
          Qwind_filtered_with, Qpv_filtered_with, Qgrid_filtered_with, Qload_filtered_with, Qsupply_filtered_with, Qin_trans_filtered_with, Qout_trans_filtered_with, Qloss_filtered_with, Qloss_LVline_filtered_with, Qloss_MVline_filtered_with, Qloss_trans_filtered_with,...
          Iwind_filtered_with, Ipv_filtered_with, Igrid_filtered_with, Iload_filtered_with, Isupply_filtered_with, Iin_trans_rms_filtered_with, Iout_trans_rms_filtered_with, ...
          f_wind_filtered_with, f_pv_filtered_with, f_grid_filtered_with, f_load_filtered_with, f_supply_filtered_with, ...
          Vwind_filtered_with, Vpv_filtered_with, Vgrid_filtered_with, Vload_filtered_with, V_supply_filtered_with, Vin_trans_filtered_with, Vout_trans_filtered_with, ...
          time_filterd_with] ...
           = filter_signals(with_demand_res)

    % If 'with_demand_res' is witht provided, retrieve it from the workspace
    if nargin < 1
        if evalin('base', 'exist(''with_demand_res'', ''var'')')
            with_demand_res = evalin('base', 'with_demand_res');
        else
            error('The variable "with_demand_res" is missing. Run the Simulink simulation first.');
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
        var_name_with = [variable_names{i}, '_with'];
        eval([var_name_with ' = with_demand_res.get(''' variable_names{i} ''');']);
    end

    time_with = with_demand_res.get('time');

    % Ensure all variables are column vectors (24×1)
    for i = 1:length(variable_names)
        var_name_with = [variable_names{i}, '_with'];
        eval([var_name_with ' = reshape(' var_name_with ', [], 1);']);
    end

    time_with = reshape(time_with, [], 1);

    % Initialize filtered outputs as column vectors (24×1)
    num_hours = 24;
    time_filtered_with = zeros(num_hours, 1);
    for i = 1:length(variable_names)
        var_name_filtered_with = [variable_names{i}, '_filtered_with'];
        eval([var_name_filtered_with ' = zeros(num_hours, 1);']);
    end

    % Process each hour separately (1–24)
    for i = 1:num_hours  
        % Find indices where time step is within the correct range
        valid_indices = (mod(time_with, 1) >= 0.85) & (mod(time_with, 1) < 0.99) & (time_with >= i-1) & (time_with < i);

        % Compute averages for all extracted variables
        for j = 1:length(variable_names)
            var_name_with = [variable_names{j}, '_with'];
            var_name_filtered_with = [variable_names{j}, '_filtered_with'];
            eval([var_name_filtered_with '(i, 1) = mean(' var_name_with '(valid_indices), ''omitnan'');']);
        end
        if any(valid_indices)
        time_filtered_with(i, 1) =round(time_with(find(valid_indices, 1, 'last')));
    else
        time_filtered_with(i, 1) = NaN;  % Handle cases where with valid time is found
    end
    end

    % Assign both raw and filtered variables to the workspace as **column vectors**
    for i = 1:length(variable_names)
        var_name_with = [variable_names{i}, '_with'];
        var_name_filtered_with = [variable_names{i}, '_filtered_with'];
        assignin('base', var_name_with, eval(var_name_with));  % Raw data
        assignin('base', var_name_filtered_with, eval(var_name_filtered_with));  % Filtered data
        assignin('base', 'time_filtered_with', time_filtered_with);
    end

    assignin('base', 'time_with', time_with);

    disp('Filtered values with "_with" and "_filtered_with" suffixes have been saved to the workspace as column vectors.');
    disp('Filtered time (time_filtered_with) saved to the workspace.');
end
