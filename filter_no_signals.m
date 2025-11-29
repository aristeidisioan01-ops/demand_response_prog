function [Pwind_filtered_no, Ppv_filtered_no, Pgrid_filtered_no, Pload_filtered_no, Psupply_filtered_no, Pin_trans_filtered_no, Pout_trans_filtered_no, Ploss_filtered_no, Ploss_LVline_filtered_no, Ploss_MVline_filtered_no, Ploss_trans_filtered_no, ...
          Qwind_filtered_no, Qpv_filtered_no, Qgrid_filtered_no, Qload_filtered_no, Qsupply_filtered_no, Qin_trans_filtered_no, Qout_trans_filtered_no, Qloss_filtered_no, Qloss_LVline_filtered_no, Qloss_MVline_filtered_no, Qloss_trans_filtered_no,...
          Iwind_filtered_no, Ipv_filtered_no, Igrid_filtered_no, Iload_filtered_no, Isupply_filtered_no, Iin_trans_rms_filtered_no, Iout_trans_rms_filtered_no, ...
          f_wind_filtered_no, f_pv_filtered_no, f_grid_filtered_no, f_load_filtered_no, f_supply_filtered_no, ...
          Vwind_filtered_no, Vpv_filtered_no, Vgrid_filtered_no, Vload_filtered_no, V_supply_filtered_no, Vin_trans_filtered_no, Vout_trans_filtered_no, ...
          time_filterd_no] ...
           = filter_signals(no_demand_res)

    % If 'no_demand_res' is not provided, retrieve it from the workspace
    if nargin < 1
        if evalin('base', 'exist(''no_demand_res'', ''var'')')
            no_demand_res = evalin('base', 'no_demand_res');
        else
            error('The variable "no_demand_res" is missing. Run the Simulink simulation first.');
        end
    end

    % Extract variables from Simulink output structure
    variable_names = {'Pwind', 'Ppv', 'Pgrid', 'Pload', 'Psupply', 'Pin_trans', 'Pout_trans', 'Ploss', 'Ploss_LVline', 'Ploss_MVline', 'Ploss_trans', ...
                      'Qwind', 'Qpv', 'Qgrid', 'Qload', 'Qsupply', 'Qin_trans', 'Qout_trans', 'Qloss', 'Qloss_LVline', 'Qloss_MVline', 'Qloss_trans', ...
                      'Iwind', 'Ipv', 'Igrid', 'Iload', 'Isupply', 'Iin_trans_rms', 'Iout_trans_rms' ...
                      'f_wind', 'f_pv', 'f_grid', 'f_load', 'f_supply', ...
                      'Vwind', 'Vpv', 'Vgrid', 'Vload', 'V_supply', 'Vin_trans', 'Vout_trans'};
                  
    % Create "_no" versions of all variables
    for i = 1:length(variable_names)
        var_name_no = [variable_names{i}, '_no'];
        eval([var_name_no ' = no_demand_res.get(''' variable_names{i} ''');']);
    end

    time_no = no_demand_res.get('time');

    % Ensure all variables are column vectors (24×1)
    for i = 1:length(variable_names)
        var_name_no = [variable_names{i}, '_no'];
        eval([var_name_no ' = reshape(' var_name_no ', [], 1);']);
    end

    time_no = reshape(time_no, [], 1);

    % Initialize filtered outputs as column vectors (24×1)
    num_hours = 24;
    time_filtered_no = zeros(num_hours, 1);
    for i = 1:length(variable_names)
        var_name_filtered_no = [variable_names{i}, '_filtered_no'];
        eval([var_name_filtered_no ' = zeros(num_hours, 1);']);
    end

    % Process each hour separately (1–24)
    for i = 1:num_hours  
        % Find indices where time step is noin the correct range
        valid_indices = (mod(time_no, 1) >= 0.85) & (mod(time_no, 1) < 0.99) & (time_no >= i-1) & (time_no < i);

        % Compute averages for all extracted variables
        for j = 1:length(variable_names)
            var_name_no = [variable_names{j}, '_no'];
            var_name_filtered_no = [variable_names{j}, '_filtered_no'];
            eval([var_name_filtered_no '(i, 1) = mean(' var_name_no '(valid_indices), ''omitnan'');']);
        end
        if any(valid_indices)
        time_filtered_no(i, 1) =round(time_no(find(valid_indices, 1, 'last')));
    else
        time_filtered_no(i, 1) = NaN;  % Handle cases where no valid time is found
    end
    end

    % Assign both raw and filtered variables to the workspace as **column vectors**
    for i = 1:length(variable_names)
        var_name_no = [variable_names{i}, '_no'];
        var_name_filtered_no = [variable_names{i}, '_filtered_no'];
        assignin('base', var_name_no, eval(var_name_no));  % Raw data
        assignin('base', var_name_filtered_no, eval(var_name_filtered_no));  % Filtered data
        assignin('base', 'time_filtered_no', time_filtered_no);
    end

    assignin('base', 'time_no', time_no);

    disp('Filtered values no "_no" and "_filtered_no" suffixes have been saved to the workspace as column vectors.');
    disp('Filtered time (time_filtered_no) saved to the workspace.');
end
