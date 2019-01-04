classdef misc_wrf_lifetime_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property-like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function workdir = workspace_dir
            mydir = fileparts(mfilename('fullpath'));
            workdir = fullfile(mydir, 'Workspaces', 'WRFData');
        end
        
        function filename = wrf_avg_name(years_in, variable)
            years_str = strjoin(sprintfmulti('%d', years_in), '_');
            filename = sprintf('Summer_WRF_avg_%s_%s.mat', years_str, variable);
            filename = fullfile(misc_emissions_analysis.emis_wrf_dir, filename);
        end
        
        function savename = diag_lifetime_savename(start_dates, end_dates)
            dvec = make_datevec(start_dates, end_dates);
            savename = sprintf('WRF_diagnostic_lifetime_%s_to_%s.mat', datestr(dvec(1),'yyyy-mm-dd'), datestr(dvec(end), 'yyyy-mm-dd'));
            savename = fullfile(misc_wrf_lifetime_analysis.workspace_dir, savename);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generative Methods %
        %%%%%%%%%%%%%%%%%%%%%%
        function make_wrf_averages(avg_year, varargin)
            E = JLLErrors;
            
            p = advInputParser;
            p.addParameter('variables',{'pres','no','no2','ho','vocr'});
            p.parse(varargin{:});
            pout = p.Results;
            variables = pout.variables;
            
            processing = struct();
            if ismember('vocr', variables)
                processing.vocr = misc_wrf_lifetime_analysis.setup_vocr_calc(avg_year);
            end
            if ismember('phox', variables)
                processing.phox = misc_wrf_lifetime_analysis.setup_phox_calc(avg_year);
            end
            if ismember('alpha', variables)
                processing.alpha = misc_wrf_lifetime_analysis.setup_alpha_calc(avg_year);
            end
            
            start_dates = cell(size(avg_year));
            end_dates = cell(size(avg_year));
            for i_yr = 1:numel(avg_year)
                start_dates{i_yr} = datenum(avg_year(i_yr), 4, 1);
                end_dates{i_yr} = datenum(avg_year(i_yr), 9, 30);
            end
            
            averages = wrf_time_average(start_dates, end_dates, variables, 'processing', processing);
            xlon = averages.XLONG;
            xlat = averages.XLAT;
            for i_var = 1:numel(variables)
                this_var = variables{i_var};
                profiles = averages.(this_var);
                save_name = misc_wrf_lifetime_analysis.wrf_avg_name(avg_year, this_var);
                save(save_name, 'xlon', 'xlat', 'profiles');
            end
        end 
        
        function wi = load_test_wrf_tile(avg_year)
            wrf_file = find_wrf_path('us','daily',sprintf('%04d-01-01', avg_year),'fullpath');
            try 
                wi = ncinfo(wrf_file);
            catch err
                if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:unableToOpenFileforRead')
                    wrf_file = strrep(wrf_file, 'wrfout', 'wrfout_subset');
                    wi = ncinfo(wrf_file);
                else
                    rethrow(err)
                end
            end
        end
        
        function vocr_processing = setup_vocr_calc(avg_year)
            % SETUP_VOCR_CALC(AVG_YEAR) Set up the processing structure to
            % calculate VOCR in WRF_TIME_AVERAGE. Requires the year being
            % averaged to figure out what VOCs are available.
            
            % First find all reactions that involve a VOC + OH.
            % Fortunately, the R2SMH a.k.a. RACM2_Berkeley2 mechanism
            % produces "OHVOC" as a diagnostic for these reactions. 
            
            rxn_net = KPP_OOP.ReactionNetwork('r2smh');
            rxns = rxn_net.FindReactionsWithProducts('OHVOC');
            
            % Now we need to get the list of VOCs that we want to average.
            wi = misc_wrf_lifetime_analysis.load_test_wrf_file(avg_year);
            wrf_vars = {wi.Variables.Name};
            rate_const_by_voc = struct();
            missing_vocs = {};
            for i_rxn = 1:numel(rxns)
                % find the reactant other than OH
                reactants = rxns{i_rxn}.reactant_names;
                voc = reactants{~strcmpi(reactants, 'HO')};
                % is the VOC in the wrf variables (case insensitive)?
                voc_idx = strcmpi(voc, wrf_vars);
                if ~any(voc_idx)
                    missing_vocs = veccat(missing_vocs, {voc});
                else
                    rate_const_by_voc.(wrf_vars{voc_idx}) = rxns{i_rxn}.rate_handle;
                end
            end
            
            % let me know if we're missing some variables
            if ~isempty(missing_vocs)
                warning('vocr_calc:missing_vocs', '%d of %d VOCs not present in wrfout file: %s', length(missing_vocs), length(rxns), strjoin(missing_vocs, ', '));
            end
            
            % for the rate constants we also need temperature and number
            % density
            wrf_vars_needed = veccat(fieldnames(rate_const_by_voc), {'temperature', 'ndens'}, 'column');
            
            vocr_processing = struct('variables', {wrf_vars_needed}, 'proc_fxn', @calc_vocr_internal);
            
            function vocr = calc_vocr_internal(Wrf)
                fns = fieldnames(rate_const_by_voc);
                for i_fn = 1:numel(fns)
                    myvoc = fns{i_fn};
                    if i_fn == 1
                        vocr = zeros(size(Wrf.(myvoc)));
                    end
                    
                    rate_expr = rate_const_by_voc.(myvoc);
                    % assuming the concentration in WRF is ppmv, convert to
                    % number density, then compute VOCR for that species.
                    voc_nd = Wrf.(myvoc) .* 1e-6 .* Wrf.ndens;
                    vocr = vocr + voc_nd .* rate_expr(Wrf.temperature, Wrf.ndens);
                end
            end
        end
        
        function phox_processing = setup_phox_calc(avg_year)
            % Assuming that most PHOx occurs via O1D + H2O -> 2OH and HCHO
            % + hv -> 2HO2 + CO, we can calculate it by:
            %
            %   1a) Calculate O1D steady-state concentration
            %       jO1D * [O3] = ko1d*[O1D] + 2.2e-10 * [H2O] * [O1D]
            %    => [O1D] = (jO1D * [O3])/(ko1d + 2.2e-10*[H2O])
            %
            %   1b) PHOx_o3 = 2.2e-10 * [H2O] * [O1D]
            %
            %   2) Add in production from HCHO photolysis:
            %       PHOx_hcho = 2 * jCH2Or * [HCHO]
            %
            % Water concentration is not given in WRF directly. Water
            % concentrations are given in kg/kg. Assuming we have the ALT
            % variable, it gives inverse density (I hope of air!). Given
            % that, then
            %
            % [H2O] = QVAPOR (kg water/kg air) * 1/ALT (m3/kg air) *
            %   mol H2O / 18.02e-3 kg H2O * 6.022e23 molec/mol * 1 m3 /
            %   (100 cm)^3 -> molec/cm^3 H2O
            
            % get the rate constant for O1D relaxation to O3P;
            rxn_network = KPP_OOP.ReactionNetwork('r2smh');
            rxns = rxn_network.FindReactionsWithReactants('O1D');
            rxns = rxns{cellfun(@(x) isscalar(x.reactants), rxns)};
            o1d_relax_rate = rxns.rate_handle;
            h2o_rate = 2.2e-10;
            
            % todo: add photolysis const
            p_o1d = 'PHOTR_O31D';
            p_hcho = 'PHOTR_CH2OR';
            req_vars = {'ALT', 'QVAPOR', 'o3', 'hcho', p_o1d, p_hcho};
            wi = misc_wrf_lifetime_analysis.load_test_wrf_tile(avg_year);
            wrf_vars = {wi.Variables.Name};
            xx = ~ismember(req_vars, wrf_vars);
            if any(xx)
                warning('missing_var:phox_calc', 'Missing the following variables in %d for PHOx calc: %s - will not calculate PHOx', avg_year, strjoin(req_vars(xx), ', '));
                phox_processing = struct('variables', {'P'}, 'proc_fxn', @no_phox);
            else
                phox_processing = struct('variables', veccat(req_vars, {'temperature', 'ndens'}), 'proc_fxn', @calc_phox);
            end
            
            function phox = no_phox(Wrf)
                phox = nan(size(Wrf.P));
            end
            
            function phox = calc_phox(Wrf)
                density = 1 ./ Wrf.ALT;
                qvapor = Wrf.QVAPOR;
                h2o = qvapor .* density ./ 18.02e-3 .* 6.022e23 ./ 100.^3;
                
                o3 = Wrf.o3;
                j_o1d = Wrf.(p_o1d);
                j_hcho = Wrf.(p_hcho);
                k_o1d = o1d_relax_rate(Wrf.temperature, Wrf.ndens);
                
                o1d = (j_o1d .* o3) ./ (h2o_rate .* h2o + k_o1d );
                
                hcho = Wrf.hcho;
                phox = 2 .* hcho .* j_hcho + o1d .* h2o .* h2o_rate;
            end
        end
        
        function alpha_processing = setup_alpha_calc(avg_year)
            % Rather than try to calculate the steady state RO2 values, I'm
            % going to take the approach used in Perring et al. 2010 (ACP p
            % 7215) where she identifies that P(O3)/P(ANs) = (2-2a)/a,
            % assuming that each RO2 + NO would produce 2 O3 if alpha were
            % 0 (one from RO2 + NO and the second from HO2 + NO).
            % Fortunately, for most years I have LNOXA (i.e. P(ANs)) and
            % PO3 saved in the WRF output!
            
            req_vars = {'LNOXA', 'PO3'};
            wi = misc_wrf_lifetime_analysis.load_test_wrf_tile(avg_year);
            wrf_vars = {wi.Variables.Name};
            xx = ~ismember(req_vars, wrf_vars);
            if any(xx)
                warning('missing_var:alpha_calc', 'Missing the following variables in %d for PHOx calc - will not calculated PHOx', avg_year, strjoin(req_vars(xx), ', '));
                alpha_processing = struct('variables', {'P'}, 'proc_fxn', @no_alpha);
            else
                alpha_processing = struct('variables', req_vars, 'proc_fxn', @calc_alpha);
            end
            
            function alpha = no_alpha(Wrf)
                alpha = nan(size(Wrf.P));
            end
            
            function alpha = calc_alpha(Wrf)
                alpha = Wrf.LNOXA ./ Wrf.PO3;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Utility Methods %
        %%%%%%%%%%%%%%%%%%%
        
        function wi = load_test_wrf_file(wrf_year)
            wrf_file = find_wrf_path('us','daily',sprintf('%04d-01-01', wrf_year),'fullpath');
            try
                wi = ncinfo(wrf_file);
            catch err
                if strcmpi(err.identifier, '?')
                    wrf_file = strrep(wrf_file, 'wrfout', 'wrfout_subset');
                    wi = ncinfo(wrf_file);
                else
                    rethrow(err)
                end
            end
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
        
        function [profiles, xlon, xlat] = load_wrf_profiles_for_years(years, specie, varargin)
            p = advInputParser;
            p.addParameter('avg_levels', []);
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_levels = pout.avg_levels;
            
            avg = RunningAverage();
            for i_yr = 1:numel(years)
                Profs = load(misc_wrf_lifetime_analysis.wrf_avg_name(years(i_yr), specie));
                if i_yr == 1
                    xlon = double(Profs.xlon);
                    xlat = double(Profs.xlat);
                end
                
                if ~isempty(avg_levels)
                    prof_data = nanmean(Profs.profiles(:,:,avg_levels), 3);
                else
                    prof_data = Profs.profiles;
                end
                
                % unlike the BEHR VCDs in the main emissions class, the WRF
                % data isn't weighted by area (since all grid cells should
                % have about the same area) so to go from 1 year to
                % multi-year averages, we can just do a simple
                % average-of-averages.
                avg.addData(double(prof_data));
            end
            profiles = avg.getWeightedAverage();
        end
        
        function [vcds, lon, lat] = compute_wrf_vcds_for_years(years, specie)
            E = JLLErrors;
            [profs, lon, lat] = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(years, specie);
            profs = profs * 1e-6; % assume profiles are in ppm
            [pres, plon, plat] = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(years, 'pres');
            
            if ~isequal(plon, lon) || ~isequal(plat,lat)
                E.callError('grid_mismatch', 'Profiles and pressure defined on different grids')
            end
            
            % Put the third dimension first, so that iterating over the
            % second and third dims goes over all profiles.
            if ndims(profs) ~= 3 || ndims(pres) ~= 3
                E.notimplemented('profiles or pressures are not 3D')
            end
            
            profs = shiftdim(profs, 2);
            pres = shiftdim(pres, 2);
            
            sz = size(profs);
            vcds = nan(sz(2:end));
            n_vcd = numel(vcds);
            wb = waitbar(0,'Computing WRF VCDs');
            for i_prof = 1:n_vcd
                waitbar(i_prof/n_vcd, wb);
                vcds(i_prof) = integPr2(profs(:,i_prof), pres(:,i_prof), pres(1,i_prof));
            end
            delete(wb);
        end
        
        function [wrf_lon, wrf_lat, wrf_data] = wrf_data_in_box(center_lon, center_lat, box_radius, wrf_lon, wrf_lat, wrf_data)
            xx_old = [];
            yy_old = [];
            
            [x_center, y_center] = find_closest_ind(center_lon, center_lat);
            xx = abs(wrf_lon(:, y_center) - center_lon) < box_radius;
            yy = abs(wrf_lat(x_center,:) - center_lat) < box_radius;
            
            while ~isequal(xx, xx_old) || ~isequal(yy, yy_old)
                % Since the WRF data isn't in an equirectangular grid, the
                % x-indices that give the right width at the top of the box
                % might not at the bottom, and likewise for the y-indices.
                % What this loop does is iterate until it finds indices
                % that give at least the requested width at all edges of
                % the box. It has to iterate since the x-indices chosen
                % will affect the width spanned by the y-indices and
                % vice-versa.
                xx_old = xx;
                yy_old = yy;
                
                xx1 = abs(wrf_lon(:, find(yy,1,'first')) - center_lon) < box_radius;
                xx2 = abs(wrf_lon(:, find(yy,1,'last')) - center_lon) < box_radius;
                yy1 = abs(wrf_lat(find(xx,1,'first'), :) - center_lat) < box_radius;
                yy2 = abs(wrf_lat(find(xx,1,'last'), :) - center_lat) < box_radius;
                
                xx = xx1 | xx2;
                yy = yy1 | yy2;

            end
            
            wrf_lon = wrf_lon(xx,yy);
            wrf_lat = wrf_lat(xx,yy);
            wrf_data = wrf_data(xx,yy,:);
            
            function [xc, yc] = find_closest_ind(target_lon, target_lat)
                r = (wrf_lon(:) - target_lon(:)).^2 + (wrf_lat(:) - target_lat).^2;
                [~, imin] = min(r);
                [xc, yc] = ind2sub(size(wrf_lon), imin);
            end
        end
        
        function locs = convert_locs(locs)
            E = JLLErrors;
            if iscellstr(locs)
                loc_inds = misc_emissions_analysis.convert_input_loc_inds(locs);
                locs = misc_emissions_analysis.read_locs_file();
                locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            elseif isnumeric(locs)
                loc_inds = locs;
                locs = misc_emissions_analysis.read_locs_file();
                locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            elseif ~isstruct(locs)
                E.badinput('LOCS must be a cell array of strings/chars, a numeric array, or a structure')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % Plotting Methods %
        %%%%%%%%%%%%%%%%%%%%
        function [oh_conc, oh_std, radii] = sample_wrf_conc_radii_by_year(locs, years, specie)
            E = JLLErrors;
            locs = misc_wrf_lifetime_analysis.convert_locs(locs);
            
            radii = 0:0.1:2;
            angles = linspace(0, 2*pi, 72);
            
            n_years = numel(years);
            n_locs = numel(locs);
            n_radii = numel(radii);
            
            oh_conc = nan(n_radii, n_years, n_locs);
            oh_std = nan(n_radii, n_years, n_locs);
            
            for i_yr = 1:n_years
                y = years(i_yr);
                window = (y-1):(y+1);
                [WRF.oh, WRF.lon, WRF.lat] = misc_wrf_lifetime_analysis.load_wrf_profiles_for_years(window, specie, 'avg_levels', 1:5);
                WRF.oh = WRF.oh * 1e-6; % convert ppm to pure mixing ratio
                for i_loc = 1:n_locs
                    max_radius = unique(locs(i_loc).BoxSize(3:4));
                    if numel(max_radius) > 1
                        E.notimplemented('Did not expect the box width to be different on opposite sides')
                    end
                    
                    % What we want is to compute the average OH at
                    % different distances from the city center, so we'll
                    % create an interpolate for a small box around the city
                    % and sample at a number of points in circles with
                    % different radii for the average.
                    this_loc = locs(i_loc);
                    [loc_lon, loc_lat, loc_oh] = misc_wrf_lifetime_analysis.wrf_data_in_box(this_loc.Longitude, this_loc.Latitude, max_radius*2,...
                        WRF.lon, WRF.lat, WRF.oh);
                    loc_interp = scatteredInterpolant(loc_lon(:), loc_lat(:), loc_oh(:));
                    
                    for i_r = 1:n_radii
                        radius = radii(i_r);
                        if radius > max_radius
                            continue
                        end
                        if radius == 0
                            oh_conc(i_r, i_yr, i_loc) = loc_interp(this_loc.Longitude, this_loc.Latitude);
                        else
                            sample_lons = this_loc.Longitude + radius .* cos(angles);
                            sample_lats = this_loc.Latitude + radius .* sin(angles);
                            
                            sample_oh = loc_interp(sample_lons, sample_lats);
                            oh_conc(i_r, i_yr, i_loc) = nanmean(sample_oh);
                            oh_std(i_r, i_yr, i_loc) = nanstd(sample_oh);
                        end
                    end
                end
            end
        end
        
        function figs = plot_wrf_conc_radii(locs, specie)
            
            locs = misc_wrf_lifetime_analysis.convert_locs(locs);
            years = 2006:2013;
            [oh_conc, oh_stds, radii] = misc_wrf_lifetime_analysis.sample_wrf_conc_radii_by_year(locs, years, specie);
            
            n_locs = numel(locs);
            n_years = numel(years);
            line_cols = nan(n_years,3);
            cmap = 'jet';
            for i_yr = 1:n_years
                line_cols(i_yr,:) = map2colmap(i_yr, 1, n_years, cmap);
            end
            
            figs = gobjects(n_locs,1);
            for i_loc = 1:n_locs
                this_loc = locs(i_loc);
                loc_lat = this_loc.Latitude;
                % convert radii to distance in kilometers
                distance = nan(size(radii));
                for i_r = 1:numel(radii)
                    distance(i_r) = m_idist(0, loc_lat, radii(i_r), loc_lat)/1000;
                end
                figs(i_loc) = figure;
                %l = gobjects(n_years,1);
                for i_yr = 1:n_years
                    this_oh_conc = oh_conc(:,i_yr,i_loc)*2e19; % convert mixing ratio to approx number density
                    this_oh_std = oh_stds(:,i_yr,i_loc)*2e19;
                    line(distance, this_oh_conc, 'linewidth', 2, 'marker', 'o', 'color', line_cols(i_yr,:),...
                        'markerfacecolor', line_cols(i_yr,:), 'markeredgecolor', 'k');
                    scatter_errorbars(distance, this_oh_conc, this_oh_std, 'color', line_cols(i_yr,:));
                end
                
                rr_has_data = any(~isnan(oh_conc(:,:,i_loc)),2);
                l1 = find(rr_has_data,1,'first');
                l2 = find(rr_has_data,1,'last')+2;
                lim_radii = veccat(distance(1) - mean(diff(distance)), distance, distance(2) + mean(diff(distance)));
                xlim([lim_radii(l1), lim_radii(l2)]);
                
                cb = colorbar;
                caxis([min(years), max(years)]);
                colormap(cmap);
                title(this_loc.Location);
                xlabel('Distance from city center (km)');
                ylabel(sprintf('[%s] (molec. cm^{-3})', upper(specie)));
                set(gca,'fontsize',16)
            end
        end
    end
end

