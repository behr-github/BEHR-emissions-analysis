classdef misc_wrf_lifetime_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        
        function workdir = workspace_dir
            mydir = fileparts(mfilename('fullpath'));
            workdir = fullfile(mydir, 'Workspaces', 'WRFData');
        end
        
        function savename = diag_lifetime_savename(start_dates, end_dates)
            dvec = make_datevec(start_dates, end_dates);
            savename = sprintf('WRF_diagnostic_lifetime_%s_to_%s.mat', datestr(dvec(1),'yyyy-mm-dd'), datestr(dvec(end), 'yyyy-mm-dd'));
            savename = fullfile(misc_wrf_lifetime_analysis.workspace_dir, savename);
        end
        
        function average_wrf_diag_lifetime(start_dates, end_dates, varargin)
            p = advInputParser;
            p.addParameter('hours', 13:23);
            p.addParameter('overwrite', nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            wrf_hours = pout.hours;
            do_overwrite = pout.overwrite;
            
            savename = misc_wrf_lifetime_analysis.diag_lifetime_savename(start_dates, end_dates);
            if exist(savename, 'file')
                if ~opt_ask_yn(sprintf('%s exists, overwrite?', savename), do_overwrite, '"overwrite"')
                    fprintf('%s exists, aborting\n', savename);
                    return
                end
            end
            
            dvec = make_datevec(start_dates, end_dates);
            
            tau_direct = [];
            tau_inverse = [];
            nox = [];
            loss = [];
            
            for i_hour = 1:numel(wrf_hours)
                this_hr = wrf_hours(i_hour);
                avg_tau_direct = RunningAverage();
                avg_tau_inverse = RunningAverage();
                avg_nox = RunningAverage();
                avg_loss = RunningAverage();
                for dnum = dvec
                    fprintf('Averaging hour %d (%s)\n', wrf_hours(i_hour), datestr(dnum));
                    wrf_file = find_wrf_path('us', 'daily', datenum(year(dnum), month(dnum), day(dnum), this_hr, 0, 0), 'fullpath');
                    [this_tau, lon, lat, this_nox, this_loss] = misc_wrf_lifetime_analysis.calculate_wrf_diag_lifetime(wrf_file, 'lifetime', 'total', 'effective', true);
                    
                    avg_tau_direct.addData(this_tau);
                    avg_tau_inverse.addData(1 ./ this_tau);
                    avg_nox.addData(this_nox);
                    avg_loss.addData(this_loss);
                end
                
                if isempty(tau_direct)
                    output_sz = [size(this_tau), numel(wrf_hours)];
                    tau_direct = nan(output_sz);
                    tau_inverse = nan(output_sz);
                    nox = nan(output_sz);
                    loss = nan(output_sz);
                end
                
                if ndims(this_tau) ~= 3
                    error('expected this_tau to have 3 dimensions')
                end
                
                % Alright, so for calculating average lifetime, there's two
                % arguments I can make for which way to do the averaging:
                %
                %   1) the line densities represent a normal, statistical
                %   average of the exponential decays, so just do a normal
                %   average of the lifetimes. But, that isn't quite right,
                %   because the average of the exponentials isn't equal to
                %   an exponential with an average lifetime:
                %
                %   2) lifetimes add inversely, when you have parallel loss
                %   routes, so we should add the lifetimes inversely in the
                %   average as well. Mathematically, for the property that
                %   the average of n identical values, v, must be equal to
                %   v to work out, this requires taking the average of the
                %   inverses then inverting that. This one isn't a perfect
                %   answer either though, because the loss processes aren't
                %   really happening in parallel.
                %
                %   3) Keep track of the average [NO_x] and loss rate, then
                %   divide them at the end to get the "lifetime of the
                %   average". 
                %
                % I'll save the data for all three since this takes a while
                % to run.
                
                tau_direct(:, :, :, i_hour) = avg_tau_direct.getWeightedAverage();
                tau_inverse(:, :, :, i_hour) = 1 ./ avg_tau_inverse.getWeightedAverage();
                nox(:, :, :, i_hour) = avg_nox.getWeightedAverage();
                loss(:, :, :, i_hour) = avg_loss.getWeightedAverage();
            end
            
            save(savename, 'tau_direct', 'tau_inverse', 'nox', 'loss', 'lon', 'lat', 'wrf_hours', 'dvec', '-v7.3');
        end
        
        function [tau, lon, lat, nox, total_loss_rate] = calculate_wrf_diag_lifetime(wrf_file, varargin)
            E = JLLErrors;

            p = advInputParser;
            p.addParameter('lifetime', '')
            p.addParameter('effective', nan);
            p.addParameter('wrf_inds', {});
            
            p.parse(varargin{:});
            pout = p.Results;
            
            lifetime_type = opt_ask_multichoice('Which lifetime to compute?', {'total', 'ANs', 'HNO3'}, pout.lifetime, '"lifetime"', 'list', true);
            is_effective_lifetime = opt_ask_yn('Calculate effective lifetime (loss - prod)?', pout.effective, '"effective"');
            
            wrf_inds = pout.wrf_inds;
            latlon_extra_args = {};
            read_extra_args = {};
            if ~isempty(wrf_inds)
                [start, count, stride] = ncind2scs(wrf_inds{:}, 'min_ndim', 3);
                latlon_extra_args(end+1:end+3) = {start, count, stride};
                [start, count, stride] = ncind2scs(wrf_inds{:}, 'min_ndim', 4);
                read_extra_args(end+1:end+3) = {start, count, stride};
            end
            
            switch lower(lifetime_type)
                case 'total'
                    loss_diag = {'LNOXA','LNOXHNO3','LNOXPAN','LNOXB','LNOX'};
                    prod_diag = {'PNOXD', 'PNOXON', 'PNOXHNO3', 'PNOXPAN'};
                case 'ans'
                    loss_diag = {'LNOXA'};
                    prod_diag = {'PNOXD', 'PNOXON'}; % PNOXD is recycling of isoprene peroxides; PNOXON is phololysis of organic nitrates
                case 'hno3'
                    loss_diag = {'LNOXHNO3'};
                    prod_diag = {'PNOXHNO3'};
                otherwise 
                    E.notimplemented('No production/loss diagnostics defined for LIFETIME_TYPE = "%s"', lifetime_type)
            end
            
            lon = ncread(wrf_file, 'XLONG', latlon_extra_args{:});
            lat = ncread(wrf_file, 'XLAT', latlon_extra_args{:});
            nox = ncread(wrf_file, 'no', read_extra_args{:}) + ncread(wrf_file, 'no2', read_extra_args{:});
            
            total_loss_rate = zeros(size(nox));
            for i_loss = 1:numel(loss_diag)
                total_loss_rate = total_loss_rate + ncread(wrf_file, loss_diag{i_loss}, read_extra_args{:});
            end
            if is_effective_lifetime
                for i_prod = 1:numel(prod_diag)
                    total_loss_rate = total_loss_rate - ncread(wrf_file, prod_diag{i_prod}, read_extra_args{:});
                end
            end
            
            % Constrain the loss rate to be >= 0. If it is <= 0, this will
            % give an infinite lifetime, which is appropriate b/c if
            % production outweighs loss, then the lifetime is effectively
            % infinite.
            total_loss_rate = max(total_loss_rate, 0);
            % According to Azimeh, the diagnostics are in ppm/s, so we can
            % divide the [NOx] and loss rate directly to get lifetime in
            % seconds (then convert to hours).
            tau = nox ./ total_loss_rate / 3600;
        end
        
        function calculate_wrf_inst_lifetime(wrf_file)
            % To calculate the instantaneous NOx lifetime, first we need to
            % find all the NOx loss and production reactions and calculate
            % the net loss rate. Since many of the reactants will be peroxy
            % radical intermediates, we will need to calculate their steady
            % state concentrations.
            
            r2smh = KPP_OOP.ReactionNetwork('r2smh');
            loss_diagnostics = {'LNOXHNO3', 'LNOXA', 'LNOXB', 'LNOXPAN', 'LNOX'};
            prod_diagnostics = {'PNOXHNO3', 'PNOX', 'PNOXD', 'PNOXON', 'PNOXPAN'};
            
            % "3d" plus the time dimension
            [start_3d, count_3d, stride_3d] = ncind2scs('all', 'all', 1, 'all');
            [start_2d, count_2d, stride_2d] = ncind2scs('all', 'all', 'all');
            
            wrf_conc = struct('NO', [], 'NO2', []);
            wrf_lon = ncread(wrf_file, 'XLONG', start_2d, count_2d, stride_2d);
            wrf_lat = ncread(wrf_file, 'XLAT', start_2d, count_2d, stride_2d);
            wrf_ndens = read_wrf_preproc(wrf_file, 'number density', start_3d, count_3d, stride_3d);
            wrf_t = read_wrf_preproc(wrf_file, 'temperature', start_3d, count_3d, stride_3d);
            
            
            
            loss_rxns = collect_reactants(loss_diagnostics);
            prod_rxns = collect_reactants(prod_diagnostics);
            load_reactants(wrf_file);
            
            total_loss_rate = calc_total_rate(loss_rxns, wrf_conc);
            total_prod_rate = calc_total_rate(prod_rxns);
            
            tau_nox = (wrf_conc.NO + wrf_conc.NO2) ./ (total_loss_rate - total_prod_rate); 
            
            function diag_rxns = collect_reactants(diagnostics)
                diag_rxns = {};
                for i_diag = 1:numel(diagnostics)
                    this_diag = diagnostics{i_diag};
                    % First, find all the reactions that result in this
                    % diagnostic. We need to make
                    these_rxns = r2smh.FindReactionsWithProducts(this_diag);
                    diag_rxns = veccat(diag_rxns, these_rxns);
                    for i_rxn = 1:numel(these_rxns)
                        reactants = these_rxns{i_rxn}.reactant_names;
                        for i_reactant = 1:numel(reactants)
                            this_reactant = reactants{i_reactant};
                            if ~isfield(wrf_conc, reactants{i_reactant})
                                wrf_conc.(this_reactant) = [];
                            end
                        end
                    end
                end
            end
            
            function load_reactants(wrf_file)
                % Load all the species that we can directly. Keep a list of
                % the rest that we'll have to do a steady-state calculation
                % for.
                ss_species = cell(1,100);
                i_ss = 0;
                conc_fns = fieldnames(wrf_conc);
                wrf_info = ncinfo(wrf_file);
                wrf_vars = {wrf_info.Variables.Name};
                for i_specie = 1:numel(conc_fns)
                    this_specie = conc_fns{i_specie};
                    if ismember(lower(this_specie), wrf_vars)
                        val = ncread(wrf_file, lower(this_specie), start_3d, count_3d, stride_3d);
                        unit = ncreadatt(wrf_file, lower(this_specie), 'units');
                        wrf_conc.(this_specie) = convert_units(val, unit, 'ppp') .* wrf_ndens;
                    elseif ismember(this_specie, wrf_vars)
                        val = ncread(wrf_file, this_specie, start_3d, count_3d, stride_3d);
                        unit = ncreadatt(wrf_file, this_specie, 'units');
                        wrf_conc.(this_specie) = convert_units(val, unit, 'ppp') .* wrf_ndens;
                    else
                        i_ss = i_ss + 1;
                        ss_species{i_ss} = this_specie;
                    end
                end
                
                % These are species I did not save, but which probably
                % cannot be considered in steady state. Other ones that I
                % am not sure of are: BAL2, CHO, MCTO which react with NO2
                % to form nitrates. EPX (an epoxide) may be able to be
                % considered in steady state.
                not_ss_species = {'H2O', 'ALD', 'UALD', 'MGLY', 'PHEN', 'CSL', 'EPX', 'MCT'};
                wrf_conc = rmfield(wrf_conc, not_ss_species);
                
                % The steady state concentrations can be returned directly
                % in molec. cm^{-3}
                ss_species = ss_species(1:i_ss);
                ss_species = ss_species(~ismember(ss_species, not_ss_species));
                ss_conc = r2smh.CalculateSpeciesSteadyState(ss_species, wrf_file, 'return_unit', 'number density', 'start_count_stride', {start_3d, count_3d, stride_3d});
                wrf_conc = copy_structure_fields(ss_conc, wrf_conc);
            end
            
            function rate = calc_total_rate(rxns, wrf_conc)
                % this will actually go through the list of reactions and
                % calculate the rate for each.
                conc_fns = fieldnames(wrf_conc);
                rate = zeros(size(wrf_conc.(conc_fns{1})));
                for i_rxn = 1:numel(rxns)
                    this_rxn = rxns{i_rxn};
                    if ~all(ismember(this_rxn.reactant_names, conc_fns))
                        fprintf('Skipping reaction b/c do not have one or more of the reactants: %s\n', char(this_rxn));
                    end
                    
                    if ~this_rxn.is_photolysis
                        this_rate = this_rxn.rate_handle(wrf_t, wrf_ndens);
                    else
                        wrf_datetime = date_from_wrf_filenames(wrf_file);
                        this_rate = KPP_OOP.interp_wrf_photolysis(this_rxn.rate_handle, wrf_datetime, wrf_lon, wrf_lat); 
                    end
                    
                    for i_react = 1:numel(this_rxn.reactant_names)
                        this_reactant = this_rxn.reactant_names{i_react};
                        this_rate = this_rate .* wrf_conc.(this_reactant);
                    end
                    
                    rate = rate + this_rate;
                end
                
            end
            
            
        end
        
    end
end

