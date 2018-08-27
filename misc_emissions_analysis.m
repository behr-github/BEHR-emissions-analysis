classdef misc_emissions_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties(Constant = true, Access = protected)
        % This is used to help check if all the required Git repos are in
        % the proper state. The first time any method calls the Git
        % verification method, if it passes, this is set to true so that it
        % doesn't need to be checked again.
        git_check_complete = false;
        
        % These define which field in the locations spreadsheet/structure
        % to use for which wind directions to reject
        wind_reject_field_std = 'WindRejects';
        wind_reject_field_wrf = 'WRFWindRejects';
        
        allowed_fit_types = {'lu','convolution'};
        
        % This is the standard fast/slow separation (in meters/second) used
        % if generating slow and fast line densities for the convolution
        % approach.
        fast_slow_sep = 3;
        
        dow_markers = struct('UMTWRFS', struct('marker', 'o', 'name', 'All days', 'used', false),...
            'TWRF', struct('marker', '^', 'name', 'Weekdays', 'used', false),...
            'US', struct('marker', 'h', 'name', 'Weekends', 'used', false));
            
        time_period_colors = struct('beg_2yr', struct('color', [0 0.5 0], 'name', '2005/07', 'used', false),...
            'beginning', struct('color', 'b', 'name', '2008*', 'used', false),...
            'end_2yr', struct('color', [0.5 0 0.5], 'name', '2012/13', 'used', false),...
            'end', struct('color', 'r', 'name', '2013*', 'used', false));
    end
    
    methods(Static = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function value = workspace_dir()
            my_dir = fileparts(mfilename('fullpath'));
            value = fullfile(my_dir, 'Workspaces');
        end
        
        function value = avg_save_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'SimpleAvgs');
        end
        
        function value = site_info_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'SiteData');
        end
        
        function value = line_density_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'LineDensities');
        end
        
        function value = emg_fit_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'EMGFits');
        end
        
        function value = emis_wrf_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'WRFData');
        end
        
        function value = table_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'Tables');
        end
        
        function value = wrf_vcd_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.emis_wrf_dir, 'WRF-VCDs');
        end
        
        function filename = avg_file_name(year_in, days_of_week, varargin)
            p = advInputParser;
            p.addOptional('species', 'NO2');
            p.parse(varargin{:});
            pout = p.Results;
            species = pout.species;
            
            years_str = strjoin(sprintfmulti('%d', year_in),'_');
            filename = sprintf('Summer_avg_%s_%s_%s.mat', species, years_str, days_of_week);
            filename = fullfile(misc_emissions_analysis.avg_save_dir, filename);
        end
        
        function filename = winds_file_name(start_date, end_date)
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('site_winds_%sto%s.mat', datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(end), 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.site_info_dir, filename);
        end
        
        function filename = line_density_file_name(start_date, end_date, by_sectors, wind_reject_filtered, wind_dir_weighted, use_wrf, winds_op, winds_cutoff, loc_inds, days_of_week)
            if by_sectors
                sectors_string = 'sectors';
            else
                sectors_string = 'rotated';
            end
            
            if wind_reject_filtered
                filtered_string = 'filtered';
            else
                filtered_string = 'unfiltered';
            end
            
            if wind_dir_weighted
                weighted_string = 'weighted';
            else
                weighted_string = 'unweighted';
            end
            
            if use_wrf
                data_string = 'WRF';
            else
                data_string = 'BEHR';
            end
            
            allowed_winds_ops = {'lt', 'gt'};
            if ~ismember(winds_op, allowed_winds_ops)
                E.badinput('WINDS_OP must be one of %s', strjoin(allowed_winds_ops, ', '));
            end
            if ~isnumeric(winds_cutoff) || ~isscalar(winds_cutoff) || winds_cutoff < 0
                E.badinput('WINDS_CUTOFF must be a positive scalar integer')
            end
            winds_string = sprintf('winds-%s%d', winds_op, winds_cutoff);
            
            if isempty(loc_inds)
                locs_string = 'locsall';
            else
                locs_string = ['locs', sprintf_ranges(loc_inds)];
            end
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('%s_%s_%s_%s_%s_%s_no2_%sto%s_%s.mat', data_string, sectors_string, filtered_string, weighted_string, winds_string, locs_string, datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(end), 'yyyy-mm-dd'), days_of_week);
            filename = fullfile(misc_emissions_analysis.line_density_dir, filename);
        end
        
        function filename = fits_file_name(start_date, end_date, using_wrf, loc_inds, days_of_week, fit_type)
            if using_wrf
                product_string = 'WRF';
            else
                product_string = 'BEHR';
            end
            
            if isempty(loc_inds)
                locs_string = 'locsall';
            else
                locs_string = ['locs', sprintf_ranges(loc_inds)];
            end
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('%s_emg_%s_fits_%s_%sto%s_%s.mat', product_string, fit_type, locs_string, datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(end), 'yyyy-mm-dd'), days_of_week);
            filename = fullfile(misc_emissions_analysis.emg_fit_dir, filename);
        end
        
        function filename = wrf_grid_area_file()
            filename = fullfile(misc_emissions_analysis.emis_wrf_dir, 'wrfgridarea_d01');
        end
        
        function filename = wrf_avg_prof_file(start_date, end_date)
            base_filename = sprintf('WRF_avg_profs_%s_to_%s.mat', datestr(start_date{1}, 'yyyy-mm-dd'), datestr(end_date{end}, 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.emis_wrf_dir, base_filename);
        end
        
        function fulldir = subdir_prep(root_dir, varargin)
            % Use this to setup sub-output directories. It will make sure
            % that the root directory exists (if not, it errors) and then
            % make the subdirectories if the don't exist
            E = JLLErrors;
            if ~exist(root_dir, 'dir')
                E.dir_dne(root_dir)
            end
            
            fulldir = fullfile(root_dir, varargin{:});
            if ~exist(fulldir, 'dir')
                mkdir(fulldir);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Utility methods %
        %%%%%%%%%%%%%%%%%%%
        
        function verify_git_state()
            if ~misc_emissions_analysis.git_check_complete
                % This requires that validate_date and list_behr_files be able
                % to handle discontinuous date ranges and list_behr_files knows
                % about the 'all' flag.
                % Also requires that the sprintf_ranges and do_keep_day_of_week 
                % functions are available.
                G = GitChecker;
                G.addReqCommits(behr_paths.behr_utils, 'aad7763');
                G.addReqCommits(behr_paths.utils, 'fc1cef0');
                % checkState() by default will error if the repositories
                % are not in the correct state.
                G.checkState();
            end
        end
        
        function locs = read_locs_file(varargin)
            locs = read_loc_spreadsheet();
            if ~isempty(varargin)
                loc_types = {locs.SiteType};
                xx = ismember(loc_types, varargin);
                locs = locs(xx);
            end
        end
        
        function winds = load_winds_file(start_date, end_date)
            % Loads a winds file for a given start and end date and inserts
            % up-to-date box size and wind direction filtering from the
            % trend_locations.xlsx sheet.
            E = JLLErrors;
            
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            
            winds_file = misc_emissions_analysis.winds_file_name(start_date, end_date);
            winds = load(winds_file);
            
            trend_locs = misc_emissions_analysis.read_locs_file();
            
            trend_shortnames = {trend_locs.ShortName};
            
            for a=1:numel(winds.locs)
                xx_loc = strcmp(winds.locs(a).ShortName, trend_shortnames);
                if sum(xx_loc) ~= 1
                    E.callError('location_not_found', 'Could not find location "%s" defined in the winds file %s but not the trend spreadsheet', winds.locs(a).ShortName, winds_file);
                end
                
                winds.locs(a).BoxSize = trend_locs(xx_loc).BoxSize;
                winds.locs(a).WindRejects = trend_locs(xx_loc).WindRejects;
            end
        end
        
        function dvec = make_datevec(start_date, end_date)
            % Create a date vector that enumerates all dates requested,
            % even if those are over non-continuous ranges.
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            
            if numel(start_date) ~= numel(end_date)
                E.badinput('START_DATE and END_DATE must have equal numbers of elements')
            end
            
            dvec = [];
            for a=1:numel(start_date)
                dvec = veccat(dvec, start_date(a):end_date(a));
            end
        end
        
        function [start_dates, end_dates, time_period, legend_id] = select_start_end_dates(time_period, varargin)
            E = JLLErrors;
            % Returns the start and end dates as cell arrays of datenums.
            % TIME_PERIOD may be 'beginning', 'end', 'beg_2yr', or
            % 'end_2yr' specifying the standard start/end dates or a 2-by-N
            % cell array of dates, where the first row will be used as the
            % start dates and the second row as the end dates.
            p = advInputParser;
            p.addOptional('prompt', 'Which time period to use?');
            p.parse(varargin{:});
            pout = p.Results;
            
            prompt = pout.prompt;
            
            if nargin < 1 || isempty(time_period)
                time_period = ask_multichoice(prompt, {'beginning (2007-09)', 'beg_2yr (2005,07)', 'end (2012-14)', 'end_2yr (2012-13)'}, 'list', true);
                time_period = strsplit(time_period, ' ');
                time_period = time_period{1};
            end
            
            start_month = 4;
            end_month = 9;
            
            if strcmpi(time_period, 'beginning')
                start_dates = {datenum(2007, start_month, 1), datenum(2008, start_month, 1), datenum(2009, start_month, 1)};
                end_dates = {eomdate(2007, end_month), eomdate(2008, end_month), eomdate(2009, end_month)};
                legend_id = '2008*';
            elseif strcmpi(time_period, 'end')
                start_dates = {datenum(2012, start_month, 1), datenum(2013, start_month, 1), datenum(2014, start_month, 1)};
                end_dates = {eomdate(2012, end_month), eomdate(2013, end_month), eomdate(2014, end_month)};
                legend_id = '2013*';
            elseif strcmpi(time_period, 'beg_2yr')
                start_dates = {datenum(2005, start_month, 1), datenum(2007, start_month, 1)};
                end_dates = {eomdate(2005, end_month), eomdate(2007, end_month)};
                legend_id = '2006*';
            elseif strcmpi(time_period, 'end_2yr')
                start_dates = {datenum(2012, start_month, 1), datenum(2013, start_month, 1)};
                end_dates = {eomdate(2012, end_month), eomdate(2013, end_month)};
                legend_id = '2012-13';
            elseif iscell(time_period)
                start_dates = validate_date(time_period(1,:));
                end_dates = validate_date(time_period(2,:));
            else
                E.badinput('TIME_PERIOD "%s" not recognized', time_period);
            end
        end
        
        function [dow, legend_id] = select_days_of_week(dow, varargin)
            p = advInputParser;
            p.addOptional('prompt', 'Select the days of week to use');
            p.parse(varargin{:});
            pout = p.Results;
            prompt = pout.prompt;
            dow = opt_ask_multichoice(prompt, {'UMTWRFS', 'TWRF', 'US'}, dow, 'days_of_week', 'list', true);
            
            legend_ids = struct('UMTWRFS', 'All', 'TWRF', 'Weekdays', 'US', 'Weekends');
            legend_id = legend_ids.(dow);
        end
        
        function inds = find_loc_struct_inds(locs)
            all_locs = misc_emissions_analysis.read_locs_file();
            all_locs_shortnames = {all_locs.ShortName};
            inds = nan(size(locs));
            for i_loc = 1:numel(locs)
                inds(i_loc) = find(strcmp(locs(i_loc).ShortName, all_locs_shortnames));
            end
        end
        
        function [xx, yy] = find_indicies_in_box_around_point(loc, lon, lat, radius)
            % LOC must be a scalar element of the locations structure, LON
            % and LAT must be 2D arrays of longitude and latitude
            % coordinates for an NO2 average or similar 2D field. RADIUS
            % must be a scalar number of grid cells in each direction to
            % get. If omitted, defaults to 0.
            E = JLLErrors;
            
            if ~exist('radius', 'var')
                radius = 0;
            end
             
            [xx, yy] = misc_emissions_analysis.find_lat_lon_index(loc.Longitude, loc.Latitude, lon, lat);
            
            xx = (xx - radius):(xx + radius);
            yy = (yy - radius):(yy + radius);
        end
        
        function [xx, radius] = find_indices_in_radius_around_loc(loc, lon, lat, varargin)
            % The Radius field is a carry over from Russell et al. 2012.
            % Now we use the boxes defined for the EMG fitting for
            % consistency. Alternately, you may specify your own radius, in
            % degrees
            p = advInputParser;
            p.addOptional('radius', []);
            p.parse(varargin{:});
            pout = p.Results;
            user_radius = pout.radius;
            
            if isempty(user_radius)
                radius = mean(loc.BoxSize(3:4));
            else
                radius = user_radius;
            end
            r = sqrt((lon - loc.Longitude).^2 + (lat - loc.Latitude).^2);
            xx = r < radius;
        end
        
        function [wrf_files, F] = closest_wrf_file_in_time(date_in, F)
            % Finds the WRF files closest in time to each swath in the BEHR
            % file for the DATE_IN. Returns the list of files as a cell
            % array. Can pass F in, which should be a structure returned
            % from DIRFF() of all relevant WRF files, which will speed up
            % this method.
            Data = load_behr_file(date_in, 'monthly', 'us'); % we only care about Time, which is the same in both monthly and daily products
            wrf_files = cell(size(Data));
            wrf_dir = find_wrf_path('us','daily',date_in);
            if ~exist('F','var')
                F = dirff(fullfile(wrf_dir, 'wrfout*'));
            end
            wrf_dates = date_from_wrf_filenames(F);
            for a=1:numel(Data)
                utc_datenum = omi_time_conv(nanmean(Data(a).Time(:)));
                [~, i_date] = min(abs(wrf_dates - utc_datenum));
                wrf_files{a} = F(i_date).name;
            end
        end
        
        function [lon, lat] = rotate_lon_lat(lon, lat, center_lon, center_lat, theta)
            R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
            for a=1:numel(lon)
                rot_coords = R * [lon(a) - center_lon; lat(a) - center_lat];
                lon(a) = rot_coords(1) + center_lon;
                lat(a) = rot_coords(2) + center_lat;
            end
        end
        
        function wind_logical = set_wind_conditions(location, speed_cutoff, winds_op, wind_reject_field)
            E = JLLErrors;
            if ~isstruct(location) || ~isscalar(location) || any(~isfield(location, {'ShortName', 'WindDir', 'WindSpeed'}))
                E.badinput('LOCATION must be a scalar structure with fields "ShortName", "WindDir", and "WindSpeed"')
            end
            
            if ~isnumeric(speed_cutoff) || ~isscalar(speed_cutoff) || speed_cutoff < 0
                E.badinput('SPEED_CUTOFF must be a scalar, positive number')
            end
            
            allowed_fast_slow = {'lt', 'gt'};
            if ~ismember(winds_op, allowed_fast_slow)
                E.badinput('WINDS_OP must be one of: %s', strjoin(allowed_fast_slow));
            end
            
            if ~exist('wind_reject_field','var') || isempty(wind_reject_field)
                wind_reject_field = 'WindRejects';
            elseif ~ischar(wind_reject_field)
                E.badinput('WIND_REJECT_FIELD must be a character array');
            elseif ~isfield(location, wind_reject_field) && ~strcmpi(wind_reject_field, 'none')
                E.badinput('The WIND_REJECT_FIELD "%s" is not a field in LOCATION and is not the string "none"', wind_reject_field);
            end
            
            wind_logical = true(size(location.WindDir));
            
            % Use the wind direction ranges specified in the locations
            % structure to reject wind directions with downwind
            % interferences that will cause an issue with the
            % lifetime/emissions. The ranges are an N-by-2 array, the first
            % column specifies the beginning of the range, the second
            % column the end. However, since wind directions go from -180
            % to +180, we have to handle the "wrap-around" nature of
            % angular coordinates. To do so, when the first column is less
            % than the second we use the typical "&" operation, otherwise
            % we use "|" (or) since e.g. all angles between +170 and -170
            % would be >170 or <-170.
            if ~strcmpi(wind_reject_field, 'none')
                for a=1:size(location.(wind_reject_field),1)
                    wind_dir_range = location.(wind_reject_field)(a,:);
                    if wind_dir_range(1) < wind_dir_range(2)
                        xx = location.WindDir >= wind_dir_range(1) & location.WindDir < wind_dir_range(2);
                    else
                        xx = location.WindDir >= wind_dir_range(1) | location.WindDir < wind_dir_range(2);
                    end
                    
                    wind_logical(xx) = false;
                end
            end
            
            % Handle wind speed filtering here
            if strcmpi(winds_op, 'gt')
                wind_logical(location.WindSpeed < speed_cutoff) = false;
            elseif strcmpi(winds_op, 'lt')
                wind_logical(location.WindSpeed >= speed_cutoff) = false;
            else
                E.notimplemented('WINDS_OP = %s is not implemented', winds_op);
            end
        end
        
        function emis_tau = calculate_emission_lifetime(line_dens_struct, fit_struct, wind_speed_vector)
            % First we need to compute the total uncertainty in the
            % parameters that accounts for uncertainty in the NO2 VCDs,
            % across wind integration distance, choice of wind fields, etc.
            param_uncert = calc_fit_param_uncert(fit_struct.ffit, fit_struct.param_stats.percent_ci95/100, line_dens_struct.num_valid_obs);
            
            % Then use these uncertainties to calculate the emissions and
            % lifetime and their uncertainties
            emissions_type = 'no';
            [emis_tau.emis, emis_tau.emis_uncert, emis_tau.tau, emis_tau.tau_uncert] = compute_emg_emis_tau(fit_struct.ffit.a, param_uncert(1), fit_struct.ffit.x_0, param_uncert(2), 'vec', wind_speed_vector, 'emissions_type', emissions_type);
            
            % Calculate emission and lifetime standard deviations and
            % degrees of freedom because they are needed for the two-sample
            % t-tests.
            
            param_sd = calc_fit_param_uncert(fit_struct.ffit, fit_struct.param_stats.percentsd/100, line_dens_struct.num_valid_obs);
            [~, emis_tau.emis_sd, ~, emis_tau.tau_sd] = compute_emg_emis_tau(fit_struct.ffit.a, param_sd(1), fit_struct.ffit.x_0, param_sd(2), 'vec', wind_speed_vector, 'emissions_type', emissions_type);
            
            % In the fitting, we consider the number of measurements to be
            % the number of points in the line density. Since we are
            % fitting 5 parameters, we lose 5 degrees of freedom.
            emis_tau.n_dofs = sum(~isnan(line_dens_struct.linedens)) - 5;
        end
        
        function [total_nei_no, nei_lon, nei_lat] = load_nei_by_year(nei_year)
            
            E = JLLErrors;
            
            for a=1:numel(nei_year)
                % Make the input path
                tmp_path = find_wrf_path('us','daily',datenum(nei_year(a),1,1));
                path_parts = strsplit(tmp_path, '/');
                
                % The path on my computer is something like '/Volumes/share-wrfN/...'
                % and we just want to get which network drive it should be on.
                % The first three parts of the split path should be an empty
                % string, Volumes, and the share. This will put a / at the
                % beginning
                wrf_share = strjoin(path_parts(1:4), '/');
                
                % Whichever share it's on, it should be in a consistent path
                % there - except for 2012. I was having trouble getting the
                % full year's inputs to prepare, so I had to split it into two
                % 6 month periods. The NEI emissions are the same in both, so
                % we can just pick one.
                inputs_path = fullfile(wrf_share, 'Inputs', num2str(nei_year(a)), 'IC-BC-Emis');
                if nei_year(a) == 2012
                    inputs_path = fullfile(inputs_path, 'Months01-06');
                elseif nei_year(a) == 2014
                    inputs_path = fullfile(inputs_path, 'Jan-Aug');
                end
                
                % We're going to average UTC 17-22, so we just need the second
                % 12 hr file
                nei_info = ncinfo(fullfile(inputs_path, 'wrfchemi_12z_d01'));
                input_info = ncinfo(fullfile(inputs_path, 'wrfinput_d01'));
                
                fprintf('Reading NEI data...\n');
                nei_lon = ncread(input_info.Filename, 'XLONG');
                nei_lat = ncread(input_info.Filename, 'XLAT');
                
                % Load the precalculated area - but double check that the
                % lat/lon matches the input
                area_lon = ncread(misc_emissions_analysis.wrf_grid_area_file, 'XLONG');
                area_lat = ncread(misc_emissions_analysis.wrf_grid_area_file, 'XLAT');
                if max(abs(area_lon(:) - nei_lon(:))) < 0.001 && max(abs(area_lat(:) - nei_lat(:))) < 0.001
                    grid_area = ncread(misc_emissions_analysis.wrf_grid_area_file, 'AREA');
                else
                    fprintf('Precomputed area lat/lon did not match, calculating WRF grid area...\n');
                    grid_area = wrf_grid_area(nei_lon, nei_lat);
                end
                
                nei_times = ncread(nei_info.Filename, 'Times')';
                nei_hours = hour(datenum(nei_times, 'yyyy-mm-dd_HH:MM:SS'));
                tt = nei_hours > 17 & nei_hours < 22;
                
                % Add up the emissions over the whole vertical extent; averaged
                % over 17:00 to 22:00 UTC, which is approximately the hours OMI
                % is over North America.
                nei_no = double(ncread(nei_info.Filename, 'E_NO'));
                nei_no = nansum(nanmean(nei_no(:,:,:,tt), 4), 3);
                
                % Not set up to handle NaNs 
                if any(isnan(nei_no(:))) || any(isnan(grid_area(:)))
                    E.notimplemented('Not set up to handle NaNs in NEI NO emissions or grid area');
                end
                
                % First time through the loop create the cumulative sum
                % array
                if a==1
                    total_nei_no = zeros(size(nei_no));
                end
                
                % Convert from mol NO / km^2 / hr to Mg NO / hr: molar mass of NO =
                % 30.006 g / mol = 30.006e-6 Mg / mol. Add to the running
                % sum, will normalize based on the number of years outside
                % the loop.
                total_nei_no = total_nei_no + nei_no .* grid_area .* 30.06e-6;
            end
            
            % To get the average NEI NO, we should be able to just divide
            % by the number of years that went into the calculation.
            total_nei_no = total_nei_no / numel(nei_year);
        end
        
        function fit = default_fit_structure()
            % Return a default structure skeleton for fitting information
            % to be used in cases where the fit fails but we want a
            % placeholder.
            null_value = [];
            fit.ffit = make_empty_struct_from_cell({'a','x_0','mu_x','sigma_x','B'},null_value);
            fit.emgfit = null_value;
            fit.param_stats = make_empty_struct_from_cell({'sd','percentsd','ci95','percent_ci95','r','r2'}, null_value);
            fit.f0 = null_value;
            fit.history.x = null_value;
            % Fit results is a complicated structure that I don't use, so
            % just make it an empty struct
            fit.fitresults = struct();
        end
        
        function emis_tau = default_emis_tau_structure()
            % Return a default structure skeleton for emissions and
            % lifetime information to be used in cases where the fit fails
            % but we want a placeholder.
            null_value = [];
            emis_tau = make_empty_struct_from_cell({'emis','emis_uncert','tau','tau_uncert','emis_sd','tau_sd','n_dofs','nei_emis'}, null_value);
        end
        
        function loc_inds = get_loc_inds_of_type(site_type)
            locs = misc_emissions_analysis.read_locs_file();
            loc_types = {locs.SiteType};
            
            if ~ischar(site_type)
                E.badinput('SITE_TYPE must be a character array')
            elseif strcmpi(site_type, 'all')
                loc_inds = 1:numel(locs);
                return
            elseif ~ismember(site_type, loc_types)
                E.badinput('SITE_TYPE is not a valid site type listed in the locations spreadsheet')
            end
            
            loc_inds = find(strcmp(site_type, loc_types));
        end
        
        function [loc_inds, file_loc_inds] = get_loc_inds_interactive(varargin)
            p = advInputParser;
            p.addFlag('one_loc');
            p.parse(varargin{:});
            pout = p.AdvResults;
            
            one_loc_flag = pout.one_loc;
            
            locs = misc_emissions_analysis.read_locs_file();
            loc_names = {locs.ShortName};
            loc_types = {locs.SiteType};
            loc_names = loc_names(~strcmpi(loc_types,'RuralAreas'));
            if ~one_loc_flag
                loc_inds = ask_multiselect('Choose the locations to use', loc_names, 'returnindex', true);
            else
                loc_inds = ask_multichoice('Choose the location to use', loc_names, 'returnindex', true);
            end
            
            min_ind = 1;
            max_ind = numel(loc_names);
            default_inds = 1:numel(loc_names);
            
            file_loc_inds = ask_number('Enter the location indicies in the file name', 'default', default_inds,...
                'testfxn', @(x) all(x >= min_ind & x <= max_ind), 'testmsg', sprintf('All values must be between %d and %d', min_ind, max_ind));
        end
        
        function fit_type_in = get_fit_type_interactive(varargin)
            if nargin > 0
                fit_type_in = varargin{1};
            else
                fit_type_in = '';
            end
            
            allowed_fit_types = misc_emissions_analysis.allowed_fit_types;
            if isempty(fit_type_in)
                fit_type_in = ask_multichoice('Which fitting method to use?', allowed_fit_types, 'list', true);
            elseif ~ismember(fit_type_in, allowed_fit_types)
                E.badinput('FIT_TYPE must be one of: %s', strjoin(allowed_fit_types, ', '));
            end
        end
        
        function ld_file = get_line_dens_file_interactive()
            avail_files = dir(fullfile(misc_emissions_analysis.line_density_dir, '*.mat'));
            avail_files = {avail_files.name};
            chosen_file = ask_multichoice('Choose the line density file to use', avail_files, 'list', true);
            ld_file = fullfile(misc_emissions_analysis.line_density_dir, chosen_file);
        end
        
        function [changes, loc_names, loc_coords] = collect_changes(first_time_period, second_time_period, first_weekdays, second_weekdays, varargin)
            
            p = inputParser;
            p.addParameter('loc_inds', []);
            p.addParameter('include_vcds', true);
            p.addParameter('use_wrf', false);
            p.addParameter('file_loc_inds', 1:71); % location indicies in the file name
            p.addParameter('fit_type','');
            p.parse(varargin{:});
            pout = p.Results;
            
            user_loc_inds = pout.loc_inds;
            include_vcds = pout.include_vcds;
            use_wrf = pout.use_wrf;
            file_loc_inds = pout.file_loc_inds;
            fit_type = misc_emissions_analysis.get_fit_type_interactive(pout.fit_type);
            
            [first_dates_st, first_dates_end] = misc_emissions_analysis.select_start_end_dates(first_time_period);
            [second_dates_st, second_dates_end] = misc_emissions_analysis.select_start_end_dates(second_time_period);
            
            first_file = misc_emissions_analysis.fits_file_name(first_dates_st(1), first_dates_end(end), use_wrf, file_loc_inds, first_weekdays, fit_type);
            fprintf('Loading %s\n', first_file);
            first_locs = load(first_file);
            second_file = misc_emissions_analysis.fits_file_name(second_dates_st(1), second_dates_end(end), use_wrf, file_loc_inds, second_weekdays, fit_type);
            fprintf('Loading %s\n', second_file);
            second_locs = load(second_file);
            
            first_locs.locs = misc_emissions_analysis.cutdown_locs_by_index(first_locs.locs, user_loc_inds, 'keep_order');
            second_locs.locs = misc_emissions_analysis.cutdown_locs_by_index(second_locs.locs, user_loc_inds, 'keep_order');
            loc_names = {first_locs.locs.Location};
            loc_coords.lon = [first_locs.locs.Longitude]';
            loc_coords.lat = [first_locs.locs.Latitude]';
            
            % Collect the emissions and lifetimes differences into two
            % n-by-2 arrays. Also get some metrics of the goodness of fit
            % and total mass
            
            additional_fns = {'r2','a','a_plus_B','is_fit_good'};
            emis_tau_fns = fieldnames(first_locs.locs(1).emis_tau);
            all_fns = veccat(emis_tau_fns, additional_fns, 'column');
            default_mat = nan(numel(first_locs.locs), 2);
            for i_fn = 1:numel(all_fns)
                changes.(all_fns{i_fn}) = default_mat;
                for i_loc = 1:numel(first_locs.locs)
                    if strcmpi(all_fns{i_fn}, 'a_plus_B')
                        first_value = a_plus_b(first_locs.locs(i_loc));
                        second_value = a_plus_b(second_locs.locs(i_loc));
                    elseif strcmpi(all_fns{i_fn}, 'is_fit_good')
                        first_value = misc_emissions_analysis.is_fit_good_by_loc(first_locs.locs(i_loc));
                        second_value = misc_emissions_analysis.is_fit_good_by_loc(second_locs.locs(i_loc));
                    else
                        first_value = find_substruct_field(first_locs.locs(i_loc), all_fns{i_fn});
                        second_value = find_substruct_field(second_locs.locs(i_loc), all_fns{i_fn});
                        
                        if isempty(first_value)
                            first_value = nan;
                        end
                        if isempty(second_value)
                            second_value = nan;
                        end
                    end
                    
                    changes.(all_fns{i_fn})(i_loc, :) = [first_value, second_value];
                end
            end
            
            if include_vcds
                changes.vcds = default_mat;
                changes.hcho_vcds = default_mat;
                if ~use_wrf
                    changes.vcds(:,1) = misc_emissions_analysis.avg_vcds_around_loc(first_locs.locs, first_time_period, first_weekdays);
                    changes.vcds(:,2) = misc_emissions_analysis.avg_vcds_around_loc(second_locs.locs, second_time_period, second_weekdays);
                    changes.hcho_vcds(:,1) = misc_emissions_analysis.avg_vcds_around_loc(first_locs.locs, first_time_period, first_weekdays, 'species', 'hcho');
                    changes.hcho_vcds(:,2) = misc_emissions_analysis.avg_vcds_around_loc(second_locs.locs, second_time_period, second_weekdays, 'species', 'hcho');
                else
                    changes.vcds(:,1) = misc_emissions_analysis.avg_wrf_vcds_around_loc(first_locs.locs, first_time_period, 'no2');
                    changes.vcds(:,2) = misc_emissions_analysis.avg_wrf_vcds_around_loc(second_locs.locs, second_time_period, 'no2');
                    changes.hcho_vcds(:,1) = misc_emissions_analysis.avg_wrf_vcds_around_loc(first_locs.locs, first_time_period, 'hcho');
                    changes.hcho_vcds(:,2) = misc_emissions_analysis.avg_wrf_vcds_around_loc(second_locs.locs, second_time_period, 'hcho');
                end
            end
            
            changes.Location = {first_locs.locs.Location}';
            changes.ShortName = {first_locs.locs.ShortName}';
            
            function aB = a_plus_b(locs)
                a = find_substruct_field(locs, 'a');
                B = find_substruct_field(locs, 'B');
                if isempty(a) || isempty(B)
                    aB = nan;
                else
                    aB = a + B;
                end
            end
            
        end
        
        function is_good = is_fit_good(x, linedens, fit_info, varargin)
            p = advInputParser;
            p.addParameter('any_num_pts', false);
            p.addParameter('DEBUG_LEVEL', 0);
            p.parse(varargin{:});
            pout = p.Results;
            
            allow_any_num_pts = pout.any_num_pts;
            DEBUG_LEVEL = pout.DEBUG_LEVEL;
            
            is_good = false;
            
            % Criteria 0: are there enough non-NaN points to give a good
            % fit? Testing with the box model shows that with only 31
            % points, the fitting procedure has a hard time fitting more
            % than a narrow range of lifetimes, but does better with 61.
            if (sum(~isnan(linedens)) < 60 || sum(~isnan(linedens)) > 70) && ~allow_any_num_pts
                if DEBUG_LEVEL > 0
                    fprintf('Fit rejected by too few line density points\n');
                end
                return
            end
            
            % Criteria 1: is R2 high enough
            if fit_info.param_stats.r2 < 0.8
                if DEBUG_LEVEL > 0
                    fprintf('Fit rejected by R2\n');
                end
                return
            end
            
            % Criteria 2: is there at least 1.5 lifetimes downwind of
            % the plume center?
            if misc_emissions_analysis.n_lifetimes_downwind(x, fit_info.ffit.x_0, fit_info.ffit.mu_x) < 1.5
                if DEBUG_LEVEL > 0
                    fprintf('Fit rejected by number of lifetimes downwind\n');
                end
                return
            end
            
            % Criteria 3: check for systematic bias in the fit
            if any(misc_emissions_analysis.test_fit_for_bias(linedens, fit_info.emgfit, 'window', 20))
                if DEBUG_LEVEL > 0
                    fprintf('Fit rejected by systematic bias test\n');
                end
                return
            end
            
            % Criteria 4: reject if sigma > x_0, because that suggests that
            % the width of the emission source might be obstructing the
            % lifetime
            if fit_info.ffit.sigma_x > fit_info.ffit.x_0
                if DEBUG_LEVEL > 0
                    fprintf('Fit rejected by width of emissions\n');
                end
                return
            end
            
            is_good = true;
        end
        
        function is_good = is_fit_good_by_loc(loc, varargin)
            E = JLLErrors;
            if ~isscalar(loc) || ~isstruct(loc)
                E.badinput('LOC must be a scalar structure')
            end
            
            varargin = update_params(varargin, 'DEBUG_LEVEL', 1);
            
            is_good = misc_emissions_analysis.is_fit_good(loc.no2_sectors.x, loc.no2_sectors.linedens, loc.fit_info, varargin{:});
        end
        
        function n_taus = n_lifetimes_downwind(x, x_0, mu_x)
            n_taus = (max(x) - mu_x)./x_0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Interactive utility methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function merge_linedens_files(DEBUG_LEVEL)
            E = JLLErrors;
            if ~exist('DEBUG_LEVEL', 'var')
                DEBUG_LEVEL = 2;
            end
            
            fprintf('Select the line density files to merge.\n');
            input('Press ENTER to continue','s');
            
            [files, path] = uigetfile('*.mat', 'Select files to merge', misc_emissions_analysis.line_density_dir, 'MultiSelect', 'on');
            if isequal(path,0)
                E.userCancel;
            elseif ~iscell(files)
                % files will be a character array if only one file was
                % selected; otherwise it's a cell array.
                fprintf('Must select >1 file to merge\n')
                return
            else
                % Check that all files are "sectors" or "rotated" and
                % "allwinds" or not. Also get the days of week, which
                % should be the last group of upper case characters before
                % the file extension, and some combination of U, M, T, W,
                % R, F, S.
                [is_sectors_vec, is_filtered_vec, is_weighted_vec, winds_strings, is_wrf_vec, days_of_week] = misc_emissions_analysis.extract_info_from_file_names(files);
                
                dow_check_fxn = @(d, d1) isequal(d,d1);
                
                if any(is_sectors_vec) && ~all(is_sectors_vec)
                    E.badinput('Some, but not all, of the files selected are by sectors');
                elseif any(is_filtered_vec) && ~all(is_filtered_vec)
                    E.badinput('Some, but not all, of the files selected are filtered by wind direction');
                elseif any(is_weighted_vec) && ~all(is_weighted_vec)
                    E.badinput('Some, but not all, of the files selected are weighted by wind direction counts');
                elseif any(is_wrf_vec) && ~all(is_wrf_vec)
                    E.badinput('Some, but not all, of the files selected are for BEHR data (as opposed to WRF)');
                elseif ~all(strcmp(winds_strings, winds_strings{1}))
                    E.badinput('Some, but not all, of the files selected are using all winds');
                elseif ~all(cellfun(@(x) dow_check_fxn(x, days_of_week{1}), days_of_week))
                    E.badinput('Not all of the files are for the same days of week (%s)', strjoin(days_of_week,' vs. '));
                else
                    by_sectors = all(is_sectors_vec);
                    is_filtered = all(is_filtered_vec);
                    is_weighted = all(is_weighted_vec);
                    wrf_bool = all(is_wrf_vec);
                    winds_op = regexp(winds_strings{1}, '(?<=winds\-)(lt|gt)','match','once');
                    winds_cutoff = str2double(regexp(winds_strings{1}, '(?<=winds\-[lg]t)\d','match','once'));
                    days_of_week = days_of_week{1};
                end
            end
            
            for i_file = 1:numel(files)
                if DEBUG_LEVEL > 1
                    fprintf('Loading %s\n', files{i_file});
                end
                LD(i_file) = load(fullfile(path, files{i_file}));
                
                % Check that all the date vectors are the same; we don't
                % want to merge files using different dates
                if i_file > 1 && ~isequal(LD(i_file).dvec, LD(1).dvec)
                    E.callError('datevec_mismatch', 'Different date vectors in %s and %s', files{1}, files{i_file});
                end
            end
            
            % Combine the line density structures then clear out the
            % original to save memory (can be 10+ GB).
            LD_all = make_empty_struct_from_cell(fieldnames(LD));
            % We already checked that all the date vectors are the same
            LD_all.dvec = LD(1).dvec;
            % Keep all the write dates
            LD_all.write_date = {LD.write_date};
            % Concatenate the locations, then check for and remove
            % duplicates
            LD_all.locs = veccat(LD.locs);
            clear('LD');
            
            loc_inds = misc_emissions_analysis.find_loc_struct_inds(LD_all.locs);
            
            [unique_inds, cut_down_vec] = unique(loc_inds);
            if numel(unique_inds) ~= numel(loc_inds)
                if ~ask_yn('Locations are duplicated among the files. Remove duplicates? (no will abort): ')
                    return
                end
            end
            % This will simultaneously remove duplicates and sort
            % everything.
            LD_all.locs = LD_all.locs(cut_down_vec);
            all_loc_inds = misc_emissions_analysis.find_loc_struct_inds(LD_all.locs);
            
            new_save_name = misc_emissions_analysis.line_density_file_name(LD_all.dvec(1), LD_all.dvec(end), by_sectors, is_filtered, is_weighted, wrf_bool, winds_op, winds_cutoff, all_loc_inds, days_of_week);
            if exist(new_save_name, 'file')
                if ~ask_yn(sprintf('%s exists. Overwrite? ', new_save_name))
                    return
                end
            end
            fprintf('Saving merged file as %s\n', new_save_name);
            save(new_save_name, '-v7.3', '-struct', 'LD_all');
        end
        
        function varargout = test_if_change_significant(varargin)
            p = advInputParser;
            p.addParameter('first_time_period','');
            p.addParameter('second_time_period','');
            p.addParameter('first_dow','');
            p.addParameter('second_dow','');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            [~, ~, first_time_period] = misc_emissions_analysis.select_start_end_dates(pout.first_time_period);
            first_dow = misc_emissions_analysis.select_days_of_week(pout.first_dow);
            [~, ~, second_time_period] = misc_emissions_analysis.select_start_end_dates(pout.second_time_period);
            second_dow = misc_emissions_analysis.select_days_of_week(pout.second_dow);
            
            % If we select a 2-year time period, those files only have 70
            % sites
            if regcmp(first_time_period, '2yr')
                file_inds = 1:70;
            else
                file_inds = 1:71;
            end
            if regcmp(second_time_period, '2yr') ~= regcmp(first_time_period, '2yr')
                E.notimplemented('Mix of 2 and 3 year files')
            end
            
            loc_inds = misc_emissions_analysis.get_loc_inds_interactive();
            
            [changes, loc_names] = misc_emissions_analysis.collect_changes(first_time_period, second_time_period, first_dow, second_dow, 'loc_inds', loc_inds, 'file_loc_inds', file_inds, 'include_vcds', false, 'fit_type', 'lu');
            is_significant = misc_emissions_analysis.is_change_significant(changes.tau, changes.tau_sd, changes.n_dofs);
            
            if nargout > 0
                varargout{1} = is_significant;
            else
                for i_loc = 1:numel(loc_names)
                    fprintf('%s (%s %s -> %s %s): %d\n', loc_names{i_loc}, first_time_period, first_dow, second_time_period, second_dow, is_significant(i_loc));
                end
            end
        end
        
        function loc_inds = ask_for_loc_inds(loc_inds, varargin)
            p = advInputParser;
            p.addOptional('allowed_types', {'Cities', 'PowerPlants'});
            p.parse(varargin{:});
            pout = p.Results;
            
            allowed_types = pout.allowed_types;
            
            if nargin < 1 || (isscalar(loc_inds) && isnan(loc_inds))
                ss_locs = misc_emissions_analysis.read_locs_file();
                loc_types = {ss_locs.SiteType};
                xx = ismember(loc_types, allowed_types);
                loc_names = {ss_locs(xx).Location};
                loc_inds = ask_multiselect('Select the locations to use:', loc_names, 'returnindex', true);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generation methods %
        %%%%%%%%%%%%%%%%%%%%%%
        function make_summer_averages(varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('avg_year',[]);
            p.addParameter('days_of_week','');
            p.addParameter('species','');
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_year = pout.avg_year;
            days_of_week = pout.days_of_week;
            species = pout.species;
            
            if isempty(avg_year)
                avg_year = ask_number('Enter the year (or years separated by a space) to do a summer average for', 'testfxn', @(x) all(x >= 2005 & x <= 2015), 'testmsg', 'Year(s) must be between 2005 and 2015');
            elseif ~isnumeric(avg_year) || any(avg_year < 2005 | avg_year > 2015)
                E.badinput('AVG_YEAR must be a numeric vector with values between 2005 and 2015');
            end
            
            days_of_week = misc_emissions_analysis.choose_days_of_week(days_of_week);
            species = opt_ask_multichoice('Which species to average?', {'NO2', 'HCHO'}, species, '"species"');
            
            start_date = cell(size(avg_year));
            end_date = cell(size(avg_year));
            for a=1:numel(avg_year)
                start_date{a} = datenum(avg_year(a), 4, 1);
                end_date{a} = datenum(avg_year(a), 9, 30);
            end
            
            % Make the monthly profile product average, then try to make
            % the daily one. If there's no data, it will return a NaN
            common_opts = {'DEBUG_LEVEL', 1, 'dayofweek', days_of_week};
            switch lower(species)
                case 'no2'
                    [monthly.no2, monthly.lon, monthly.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'monthly', common_opts{:});
                    [daily.no2, daily.lon, daily.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'daily', common_opts{:});
                case 'hcho'
                    [monthly.hcho, monthly.lon, monthly.lat] = omhcho_time_average(start_date, end_date);
                    daily = monthly;
            end
            
            save_name = misc_emissions_analysis.avg_file_name(avg_year, days_of_week, species);
            save(save_name, 'monthly', 'daily');
        end
        
        function make_location_winds_file(time_period, overwrite)
            % As in Laughner, Zare, and Cohen (2016, ACP) we will calculate
            % wind direction by averaging over the first 5 WRF layers in a
            % 3x3 grid centered on each location.
            
            misc_emissions_analysis.verify_git_state();
            
            if ~exist('time_period', 'var')
                time_period = '';
            end
            
            if ~exist('overwrite', 'var')
                overwrite = -1;
            end
            
            locs = misc_emissions_analysis.read_locs_file();
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
            dvec = misc_emissions_analysis.make_datevec(start_date, end_date);
            
            % Check that the save file exists
            save_file = misc_emissions_analysis.winds_file_name(start_date, end_date);
            if exist(save_file, 'file')
                % If overwrite isn't specified, ask. Otherwise, if it is
                % false and the file exists, abort.
                if overwrite < 0 && ~ask_yn(sprintf('%s already exists. Overwrite?', save_file))
                    return
                elseif ~overwrite
                    return
                end
            end
            
            wind_array = nan(numel(dvec), 6); % ndates x norbits per date (usually 4, sometimes 5, so 6 should be plenty)
            wind_cell = cell(size(wind_array));
            for a=1:numel(locs)
                locs(a).WindDir = wind_array;
                locs(a).WindSpeed = wind_array;
                locs(a).U = wind_cell;
                locs(a).V = wind_cell;
            end
            
            last_month = -1;
            last_year = -1;
            for d=1:numel(dvec)
                % We need to pick the WRF file closest in time to the OMI
                % overpass, so we will load the BEHR file for this day and
                % calculate which WRF file is closest in time to that
                % swath, then calculate wind speed and direction for each
                % file. Later, we will actually pick the swath-specific
                % wind direction for each rotation.
                
                % If we're in the same month as the last time through this
                % loop, the closest_wrf method doesn't need to get the
                % directory listing of the WRF directory again (because it
                % should be organized by month and year).
                fprintf('%s: Gathering WRF files\n', datestr(dvec(d)));
                try
                    if month(dvec(d)) == last_month && year(dvec(d)) == last_year
                        fprintf('     Using existing list of WRF files\n');
                        wrf_files = misc_emissions_analysis.closest_wrf_file_in_time(dvec(d), all_months_wrf_files);
                    else
                        fprintf('     New month: need to get the directory listing\n');
                        [wrf_files, all_months_wrf_files] = misc_emissions_analysis.closest_wrf_file_in_time(dvec(d));
                        last_month = month(dvec(d));
                        last_year = year(dvec(d));
                    end
                catch err
                    if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
                        fprintf('Cannot load file for %s, skipping\n', datestr(dvec(d)));
                        continue
                    else
                        rethrow(err)
                    end
                end
                
                % Load the bottom five layers of U and V, plus COSALPHA and
                % SINALPHA
                for a=1:numel(wrf_files)
                    fprintf('  %s, swath %d: Reading wind and lat/lon\n', datestr(dvec(d)), a);
                    U = ncread(wrf_files{a}, 'U', [1 1 1 1], [Inf, Inf, 5, Inf]);
                    V = ncread(wrf_files{a}, 'V', [1 1 1 1], [Inf, Inf, 5, Inf]);
                    cosalpha = ncread(wrf_files{a}, 'COSALPHA');
                    sinalpha = ncread(wrf_files{a}, 'SINALPHA');
                    wrf_lon = ncread(wrf_files{a}, 'XLONG');
                    wrf_lat = ncread(wrf_files{a}, 'XLAT');
                    
                    U = unstagger(U,1);
                    V = unstagger(V,2);
                    [U, V] = wrf_winds_transform(U,V,cosalpha,sinalpha);
                    
                    % As in Laughner et al. 2016, we average over the
                    % vertical layers first (c.f. misc_behr_wind_plots.m,
                    % subfunction plot_wind_magnitude_and_angle in
                    % https://github.com/behr-github/BEHR-WindEffect-analysis)
                    % Theoretically, this should be the same as averaging
                    % over all 45 values at once (I think)
                    U = nanmean(U,3);
                    V = nanmean(V,3);
                    
                    % Now loop over each location and calculate its wind
                    % for that swath
                    fprintf('  %s, swath %d: Averaging to locations\n', datestr(dvec(d)), a);
                    for l=1:numel(locs)
                        [xx,yy] = misc_emissions_analysis.find_indicies_in_box_around_point(locs(l), wrf_lon, wrf_lat, 1);
                        Ubar = nanmean(reshape(U(xx,yy),[],1));
                        Vbar = nanmean(reshape(V(xx,yy),[],1)); 
                        
                        locs(l).WindSpeed(d,a) = sqrt(Ubar.^2 + Vbar.^2);
                        locs(l).WindDir(d,a) = atan2d(Vbar,Ubar);
                        
                        % We'll also store the vertically averaged wind
                        % fields for checking later
                        locs(l).U{d,a} = U(xx,yy);
                        locs(l).V{d,a} = V(xx,yy);
                    end
                end
            end
            
            locs = misc_emissions_analysis.mark_which_winds_will_be_used(locs, dvec); %#ok<NASGU>
            
            write_date = datestr(now); %#ok<NASGU>
            
            save(save_file, 'locs', 'dvec', 'write_date');
        end
        
        function locs = mark_which_winds_will_be_used(locs, dvec)
            % Iterate through the dates used for each location and estimate
            % which winds will actually be used based on rotated_plume.
            % This is used in order to weight the slow wind speed data to
            % have the same fractional contribution to a given wind
            % direction, for which we don't just want to count up all
            % instances of wind falling in a directional bin, just those
            % that have actual satellite data.
            
            % Copied from calc_line_density on 26 Mar 2018.
            reject_details = struct('cloud_type', 'omi', 'cloud_frac', 0.2, 'row_anom_mode', 'XTrackFlags', 'check_behr_amf', true);
            
            % Initialize the WindUsedBool field for all locations
            for i_loc = 1:numel(locs)
                locs(i_loc).WindUsedBool = false(size(locs(i_loc).WindDir));
            end
            
            % Iterate over dates in the outer loop since we need to load
            % the BEHR file for each day.
            for i_date = 1:numel(dvec)
                if all(isnan(locs(1).WindSpeed(i_date,:)))
                    fprintf('No wind speeds defined for %s, so leaving all false\n', datestr(dvec(i_date)));
                    continue
                end
                fprintf('Marking useful winds for %s\n', datestr(dvec(i_date)));
                Data = load_behr_file(dvec(i_date),'daily','us');
                for i_orbit = 1:numel(Data)
                    Data(i_orbit).Areaweight = ones(size(Data(i_orbit).Longitude));
                    Data(i_orbit) = omi_pixel_reject(Data(i_orbit),'detailed',reject_details);
                    for i_loc = 1:numel(locs)
                        if strcmp(locs(i_loc).SiteType, 'RuralAreas')
                            % Rural areas do not have boxes defined, so we
                            % can't use rotate_plume here, but they aren't
                            % used in this analysis anyway.
                            continue
                        end
                        
                        xx_pixels_used = rotate_plume(Data(i_orbit), locs(i_loc).Longitude, locs(i_loc).Latitude, locs(i_loc).WindDir(i_date, i_orbit), locs(i_loc).BoxSize, 'pixels_in_box', true);
                        locs(i_loc).WindUsedBool(i_date, i_orbit) = any(Data(i_orbit).Areaweight(xx_pixels_used) > 0);
                    end
                end
            end
        end
        
        function make_rotated_line_densities(varargin)
            % MAKE_ROTATED_LINE_DENSITIES() will make the line densities
            % aligning all wind directions for all locations, asking before
            % overwriting existing files.
            %
            % MAKE_ROTATED_LINE_DENSITIES( LOC_INDICIES ) will restrict the
            % locations to those given at locs(LOC_INDICIES) in the winds
            % file (misc_emissions_analysis.winds_file_name). It will do
            % this before removing rural sites. Will still ask to overwrite
            % existing file.
            %
            % MAKE_ROTATED_LINE_DENSITIES( LOC_INDICIES, OVERWRITE ) given
            % OVERWRITE == 0, will not overwrite an existing file and
            % OVERWRITE > 0 will overwrite an existing file. OVERWRITE < 0
            % will ask before overwriting.
            if ~any(strcmp('winds_op', varargin)) && ~any(strcmp('winds_cutoff', varargin))
                varargin(end+1:end+4) = {'winds_op', 'gt', 'winds_cutoff', 3};
            end
            misc_emissions_analysis.make_line_densities(false, varargin{:});
        end
        
        function make_rotated_slow_convolution_line_densities(varargin)
            % MAKE_ROTATED_SLOW_CONVOLUTION_LINE_DENSITIES() will call
            % make_line_densities() with the following parameters enforced:
            %
            %       'weight_wind_dirs' = true
            %       'use_wind_rejects' = false
            %       'winds_op' = 'lt'
            %       'winds_cutoff' = misc_emissions_analysis.fast_slow_sep
            varargin = update_params(varargin, 'weight_wind_dirs', true, 'use_wind_rejects', false, 'winds_op', 'lt', 'winds_cutoff', misc_emissions_analysis.fast_slow_sep);
            misc_emissions_analysis.make_line_densities(false, varargin{:});
        end
        
        function make_rotated_fast_convolution_line_densities(varargin)
            % MAKE_ROTATED_FAST_CONVOLUTION_LINE_DENSITIES() will call
            % make_line_densities() with the following parameters enforced:
            %
            %       'weight_wind_dirs' = false
            %       'use_wind_rejects' = false
            %       'winds_op' = 'gt'
            %       'winds_cutoff' = misc_emissions_analysis.fast_slow_sep
            varargin = update_params(varargin, 'weight_wind_dirs', false, 'use_wind_rejects', false, 'winds_op', 'gt', 'winds_cutoff', misc_emissions_analysis.fast_slow_sep);
            misc_emissions_analysis.make_line_densities(false, varargin{:});
        end
        
        function make_sector_line_densities(varargin)
            % MAKE_SECTOR_LINE_DENSITIES() will make the line densities for
            % separate wind direction sectors for all locations, asking
            % before overwriting existing files.
            %
            % MAKE_SECTOR_LINE_DENSITIES( LOC_INDICIES ) will restrict the
            % locations to those given at locs(LOC_INDICIES) in the winds
            % file (misc_emissions_analysis.winds_file_name). It will do
            % this before removing rural sites. Will still ask to overwrite
            % existing file.
            %
            % MAKE_SECTOR_LINE_DENSITIES( LOC_INDICIES, OVERWRITE ) given
            % OVERWRITE == 0, will not overwrite an existing file and
            % OVERWRITE > 0 will overwrite an existing file. OVERWRITE < 0
            % will ask before overwriting.
            if ~any(strcmp('winds_op', varargin)) && ~any(strcmp('winds_cutoff', varargin))
                varargin(end+1:end+4) = {'winds_op', 'lt', 'winds_cutoff', 3};
            end
            misc_emissions_analysis.make_line_densities(true, varargin{:});
        end
        
        function make_line_densities(by_sectors, varargin)
            % MAKE_LINE_DENSITIES Create line density output files.
            %
            %   MAKE_LINE_DENSITIES( BY_SECTORS )
            %
            %   Not really intended to be called directly, there are other
            %   make_*_line_densities functions that will pass some of the
            %   necessary arguments automatically to make sure this is set
            %   up correctly. One required argument, BY_SECTORS, which must
            %   be a scalar logical value, indicating if line densities
            %   should be divided up into sectors (true), or rotated
            %   (false).
            %
            %   Parameters:
            %
            %   'time_period' - which time period to calculate line
            %   densities for. Must be a string recognized by
            %   select_start_end_dates(). If empty, will prompt for the
            %   time period (default).
            %
            %   'loc_indicies' - a numeric vector indicated which locations
            %   to calculate line densities for by their index in the
            %   locations spreadsheet. If empty (default) all are 
            %   calculated.
            %
            %   'do_overwrite' - scalar logical that indicates whether
            %   existing line density file should be overwritten. If -1
            %   (default), will prompt.
            %
            %   'days_of_week' - which days of week to include in the line
            %   density. Must be a string containing some subset of the
            %   characters UMTWRFS. Default is all days of week.
            %
            %   'winds_op' - the string 'gt' (default) or 'lt', which in
            %   conjunction with 'winds_cutoff' (next) determines the
            %   criteria for wind speed. 'gt' = greater than, 'lt' = less
            %   than
            %
            %   'winds_cutoff' - a scalar number indicating the wind speed
            %   criteria that goes with 'winds_op'.
            %
            %   'use_wrf' - a scalar logical indicating whether to use VCDs
            %   derived from WRF model simulation (using the
            %   preproc_AprioriVCDs utility in this repository). By
            %   default, this will use the WRFWindRejects field of the
            %   locations spreadsheet to determine which wind directions to
            %   ignore; that can be overridden by 'use_wrf_wind_rejects'.
            %   This will also not filter the line densities for cloud or
            %   row anomaly.
            %
            %   'use_wind_rejects' - a boolean that determines whether wind
            %   directions listed in the locations spreadsheet should be
            %   skipped. Default is true, false will use all wind
            %   directions.
            %
            %   'use_wrf_wind_rejects' - if rejecting by wind direction,
            %   this controls whether it uses the WindRejects (false) or
            %   WRFWindRejects (true) field in the spreadsheet. By default
            %   this is set to match 'use_wrf', but you may override that.
            %
            %   'weight_wind_dirs' - boolean, default false, that indicates
            %   whether the contribution of each wind sector to the line
            %   density should be weighted. Each sector is weighted by the
            %   number of days in it that have fast winds divided by those
            %   with slow winds; it is intended to be used to weight slow
            %   wind conditions for use with the convolution algorithm.
            
            E = JLLErrors;
            
            if ~islogical(by_sectors) || ~isscalar(by_sectors)
                E.badinput('BY_SECTORS must be a scalar logical')
            end
            
            p = advInputParser;
            p.addParameter('time_period', '');
            p.addParameter('loc_indices', []);
            p.addParameter('do_overwrite', -1);
            p.addParameter('days_of_week', 'UMTWRFS');
            p.addParameter('winds_op', 'gt')
            p.addParameter('winds_cutoff', 3);
            p.addParameter('use_wrf', false);
            p.addParameter('use_wind_rejects',true);
            p.addParameter('use_wrf_wind_rejects', nan);
            p.addParameter('weight_wind_dirs',false);
            
            p.parse(varargin{:});
            pout = p.AdvResults;
            
            misc_emissions_analysis.verify_git_state();
            
            time_period = pout.time_period;
            loc_indicies = pout.loc_indices;
            do_overwrite = pout.do_overwrite;
            days_of_week = pout.days_of_week;
            winds_op = pout.winds_op;
            winds_cutoff = pout.winds_cutoff;
            wrf_bool = pout.use_wrf;
            use_wind_rejects = pout.use_wind_rejects;
            use_wrf_rejects = pout.use_wrf_wind_rejects;
            weight_wind_dirs = pout.weight_wind_dirs;
            
            if ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('The parameter "loc_indicies" must be a numeric array with all values >= 1')
            end
            
            if (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('The parameter "do_overwrite" must be a scalar logical or number')
            end
            
            if isnan(use_wrf_rejects)
                use_wrf_rejects = wrf_bool;
            end
            
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
 
            % Find the list of BEHR files between the start and end dates
            [behr_files, behr_dir] = list_behr_files(start_date, end_date,'daily','all');
            if wrf_bool
                for i_file = 1:numel(behr_files)
                    behr_files(i_file).name = strrep(behr_files(i_file).name, 'OMI', 'WRF');
                end
                behr_dir = misc_emissions_analysis.wrf_vcd_dir;
            end
            % If we're doing line densities by sector, then we don't want
            % to reject any wind directions. We also don't want to reject
            % wind directions if explicitly told not to.
            filter_by_wind_dir = use_wind_rejects && ~by_sectors;
            if ~filter_by_wind_dir
                wind_reject_field = 'none';
            elseif use_wrf_rejects
                wind_reject_field = 'WRFWindRejects';
            else
                wind_reject_field = 'WindRejects';
            end
            
            % If overwrite not given and the save file exists, ask to
            % overwrite. Otherwise, only overwrite if instructed.
            save_name = misc_emissions_analysis.line_density_file_name(start_date, end_date, by_sectors, filter_by_wind_dir, weight_wind_dirs, wrf_bool, winds_op, winds_cutoff, loc_indicies, days_of_week);
            if exist(save_name, 'file')
                if do_overwrite < 0
                    if ~ask_yn(sprintf('%s exists. Overwrite?', save_name))
                        return
                    end
                elseif ~do_overwrite
                    return
                end
            end
           
            fprintf('Will save as %s\n', save_name);
            
            % Load the winds file with up-to-date box sizes and wind
            % sectors to reject
            winds = misc_emissions_analysis.load_winds_file(start_date, end_date);
            winds.locs = misc_emissions_analysis.append_new_spreadsheet_fields(winds.locs);
            
            
            
            % Check that the dates match up with what we're expecting (it
            % should because we load the file with those dates). Also check
            % that it matches up with the BEHR filenames.
            check_dvec = misc_emissions_analysis.make_datevec(start_date, end_date);
            if ~isequal(check_dvec, winds.dvec)
                E.callError('date_mismatch', 'Dates in winds file do not match required (%s to %s, %d dates)', datestr(start_date(1)), datestr(end_date(end)), numel(check_dvec));
            end
            
            behr_dvec = date_from_behr_filenames(behr_files);
            if ~isequal(winds.dvec(:), behr_dvec(:))
                E.callError('date_mismatch', 'Dates in the winds file (%s) do not match the BEHR files listed', winds_file)
            end

            if ~isempty(loc_indicies)
                winds.locs = misc_emissions_analysis.cutdown_locs_by_index(winds.locs, loc_indicies);
            end
            % Rural sites aren't going to be that interesting
            xx = strcmpi('Cities', {winds.locs.SiteType}) | strcmpi('PowerPlants', {winds.locs.SiteType});
            winds.locs(~xx) = [];
            % This should allow the substructure "locs" to be a sliced,
            % instead of broadcast, variable
            winds_locs_distributed = winds.locs;
            
            % Set up a default box size in case the spreadsheet has an
            % invalid box size.
            default_box = [1 2 1 1];
            
            % Should we weight the line densities by the ratio of the
            % number of times the winds go in a given direction when they
            % are slow vs. fast? If so, append the necessary arguments to
            % the call to calc_line_density. (Right now will not affect
            % sectors.)
            
            if weight_wind_dirs
                [wind_dir_weights, wind_dir_edges] = misc_emissions_analysis.calculate_wind_bin_weights(winds.locs, winds_cutoff, 'all_winds', false);
            else
                % When we hit the parfor loop, Matlab will try to transmit 
                % all variables needed in the loop to the workers. It does 
                % not care that wind_dir_weights and wind_dir_edges aren't
                % needed if weight_wind_dirs is false, so we have to give 
                % these some fill value to avoid a "variable not defined"
                % error. These are actually sliced variables in the parfor 
                % loop, so they need to be the same size as winds.locs 
                % because the parfor loop will try to send e.g. wind_dir_weights{20}
                % to the worker doing a=20.
                wind_dir_weights = cell(size(winds.locs));
                wind_dir_edges = cell(size(winds.locs));
            end
            
            for a=1:numel(winds.locs)
                opt_args = {};
                
                box_size = winds_locs_distributed(a).BoxSize;
                if any(isnan(box_size))
                    warning('NaN detected in box size. Setting to default %s for location %s', mat2str(default_box), winds_locs_distributed(a).ShortName)
                    box_size = default_box;
                end
                
                if weight_wind_dirs && ~by_sectors
                    opt_args = veccat(opt_args, {'wind_dir_weights', wind_dir_weights{a}, 'wind_weights_bins', wind_dir_edges{a}});
                end
                
                if wrf_bool
                    opt_args = veccat(opt_args, {'no_reject', true});
                end
                
                % "wind_reject_field" will have been set to 'none' if
                % either doing sectors or told explicitly not to filter by
                % wind direction. This will always filter by wind speed
                % though.
                wind_logical = misc_emissions_analysis.set_wind_conditions(winds_locs_distributed(a), winds_cutoff, winds_op, wind_reject_field);

                if by_sectors
                    fprintf('Calculating sector line densities for %s\n', winds_locs_distributed(a).ShortName);
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density_sectors(behr_dir, behr_files, winds_locs_distributed(a).Longitude, winds_locs_distributed(a).Latitude, winds_locs_distributed(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size, opt_args{:});
                else
                    fprintf('Calculating rotated line densities for %s\n', winds_locs_distributed(a).ShortName);
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density(behr_dir, behr_files, winds_locs_distributed(a).Longitude, winds_locs_distributed(a).Latitude, winds_locs_distributed(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size, 'days_of_week', days_of_week, opt_args{:});
                end
            end
            
            for a=1:numel(winds.locs)
                winds.locs(a).no2_sectors = no2(a);
            end
            
            locs = winds.locs;
            dvec = winds.dvec;
            write_date = datestr(now);
            
            save(save_name, '-v7.3', 'locs', 'dvec', 'write_date');
        end
        
        function locs = make_emg_fits(varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('time_period', '');
            p.addParameter('loc_indicies', []);
            p.addParameter('file_loc_indicies','match'); % if set to 'match' then this will be the same as loc_indicies.
            p.addParameter('add_nei', true);
            p.addParameter('days_of_week', 'UMTWRFS');
            p.addParameter('use_wrf', false);
            p.addParameter('do_overwrite', -1);
            % by default if it doesn't get it the second time, then it's
            % probably just going to randomly sample until it happens to
            % get two runs that give the same fit, which there's no reason
            % to believe that is the minimum.
            p.addParameter('max_fit_attempts', 2);  
            p.addParameter('fatal_fit_fail', true);
            p.addParameter('skip_linedens_errors', -1);
            p.addParameter('fit_type', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            time_period = pout.time_period;
            loc_indicies = pout.loc_indicies;
            file_loc_indicies = pout.file_loc_indicies;
            add_nei = pout.add_nei;
            days_of_week = pout.days_of_week;
            wrf_bool = pout.use_wrf;
            do_overwrite = pout.do_overwrite;
            max_fit_attempts = pout.max_fit_attempts;
            fatal_if_cannot_fit = pout.fatal_fit_fail;
            skip_linedens_errors = pout.skip_linedens_errors;
            fit_type_in = pout.fit_type;
            % time_period should be checked in select_start_end_dates
            
            if ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('LOC_INDICIES must be a numeric array with all values >= 1')
            end
            
            if strcmpi(file_loc_indicies, 'match')
                file_loc_indicies = loc_indicies;
            elseif ~isnumeric(file_loc_indicies) || any(file_loc_indicies(:) < 1)
                E.badinput('FILE_LOC_INDICIES must be a numeric array with all values >= 1 or the string "match"')
            end
            
            if ~isscalar(add_nei) || (~islogical(add_nei) && ~isnumeric(add_nei))
                E.badinput('ADD_NEI must be a scalar logical or numeric value');
            end
            
            if ~ischar(days_of_week)
                E.badinput('DAYS_OF_WEEK must be a character array');
            end
            
            if (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('DO_OVERWRITE must be a scalar logical or number')
            end
            
            fit_type_in = misc_emissions_analysis.get_fit_type_interactive(fit_type_in);
            
            
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
            % If overwrite not given and the save file exists, ask to
            % overwrite. Otherwise, only overwrite if instructed.
            save_name = misc_emissions_analysis.fits_file_name(start_date, end_date, wrf_bool, loc_indicies, days_of_week, fit_type_in);
            if exist(save_name, 'file')
                if do_overwrite < 0
                    if ~ask_yn(sprintf('%s exists. Overwrite?', save_name))
                        return
                    end
                elseif ~do_overwrite
                    [~,save_basename] = fileparts(save_name);
                    fprintf('%s exists already. Not overwriting\n', save_basename);
                    return
                end
            end
            
            
            if strcmpi(fit_type_in, 'convolution')
                % For the convolution approach, we want the fast line
                % densities to be unfiltered for wind direction (because in
                % theory the convolution with the slow line densities will
                % handle the downwind sources we had to filter out in the
                % normal way) and unweighted (fast line densities should
                % always be unweighted for wind diretion contribution).
                filtered_bool = false;
                weighted_bool = false;
                % However for the slow line densities we do want them
                % weighted for wind direction contribution so that they
                % match the fast line densities.
                slow_ldens_file = misc_emissions_analysis.line_density_file_name(start_date, end_date, false, filtered_bool, true, wrf_bool, 'lt', misc_emissions_analysis.fast_slow_sep, file_loc_indicies, days_of_week);
                if ~exist(slow_ldens_file, 'file')
                    [~,ldens_basename] = fileparts(slow_ldens_file);
                    % Use regular error function to have more control over the
                    % error identifier
                    error('emis_analysis:no_linedens_file', 'Slow line density file %s not found, cannot fit EMG functions', ldens_basename);
                end
                slow_line_densities = load(slow_ldens_file);
            else
                % Load the file with the line densities. For this, we never
                % want the sectors line densities (first false) and usually
                % want the file with winds greater than the separation
                % speed (3 m/s as of 27 Mar 2018).
                filtered_bool = true;
                weighted_bool = false;
                slow_line_densities = [];
            end
            ldens_file = misc_emissions_analysis.line_density_file_name(start_date, end_date, false, filtered_bool, weighted_bool, wrf_bool, 'gt', misc_emissions_analysis.fast_slow_sep, file_loc_indicies, days_of_week);
            if ~exist(ldens_file, 'file')
                [~,ldens_basename] = fileparts(ldens_file);
                % Use regular error function to have more control over the
                % error identifier
                error('emis_analysis:no_linedens_file', 'Line density file %s not found, cannot fit EMG functions', ldens_basename);
            end
            line_densities = load(ldens_file);
            
            % Check that the dates match up with what we're expecting (it
            % should because we load the file with those dates)
            check_dvec = misc_emissions_analysis.make_datevec(start_date, end_date);
            if ~isequal(check_dvec, line_densities.dvec)
                E.callError('date_mismatch', 'Dates in winds file (%s) do not match required (%s to %s)', ldens_file, datestr(start_date(1)), datestr(end_date(end)));
            end
            
            if ~isempty(loc_indicies)
                locs = misc_emissions_analysis.cutdown_locs_by_index(line_densities.locs, loc_indicies);
                if ~isempty(slow_line_densities)
                    slow_line_densities.locs = misc_emissions_analysis.cutdown_locs_by_index(slow_line_densities.locs, loc_indicies);
                end
            else
                locs = line_densities.locs;
            end
            
            % Load the NEI data. Will need to get lat/lon from
            % wrfinput_d01, b/c the wrfchemi files don't include lat-lon.
            % Would have to do some work to get this to run on the cluster.
            if add_nei
                nei_year = unique(year(check_dvec));
                [nei_avg_no, nei_lon, nei_lat] = misc_emissions_analysis.load_nei_by_year(nei_year);
            end
            % Specify even the default options so that if fit_line_density
            % changes, we know exactly what options we wanted.
            common_opts = {'fmincon_output', 'none', 'fittype', 'ssresid', 'nattempts', 20};
            if strcmpi(fit_type_in, 'convolution')
                % In Liu 2016, when she applies the convolved line
                % densities, the center offset in the exponential is set to
                % 0. I assume this is because the slow line densities
                % should already contain information about the true center
                % point of the emission.
                common_opts = veccat(common_opts, {'fixed_param','mux','fixed_val',0});
                
                % We also need to indicate that the mu_x and sigma_x
                % parameters don't matter when checking if the two fitting
                % attempts are the same
                fit_check_inds = [1 2 5];
            else
                fit_check_inds = 1:5;
            end
            
            for a=1:numel(locs)
                fprintf('Fitting %s\n', locs(a).ShortName);
                safety_count = 1;
                while true
                    if strcmpi(fit_type_in, 'convolution')
                        fit_type = convolved_fit_function(slow_line_densities.locs(a).no2_sectors.x, slow_line_densities.locs(a).no2_sectors.linedens);
                    else
                        fit_type = fit_type_in;
                    end
                    try
                        [ffit, emgfit, param_stats, f0, history, fitresults] = fit_line_density(locs(a).no2_sectors.x, locs(a).no2_sectors.linedens, 'emgtype', fit_type, common_opts{:});
                        % Try this a second time - if it gives a different
                        % answer, we should re-run, since that suggests we
                        % didn't find the minimum one time.
                        ffit_check = fit_line_density(locs(a).no2_sectors.x, locs(a).no2_sectors.linedens, 'emgtype', fit_type, common_opts{:});
                    catch err
                        msg = sprintf('Fitting %s failed with error:\n "%s"\nSkip this location and continue?', locs(a).ShortName, err.message);
                        if skip_linedens_errors > 0 || (skip_linedens_errors < 0 && ask_yn(msg))
                            locs(a).fit_info = misc_emissions_analysis.default_fit_structure;
                            break
                        else
                            rethrow(err)
                        end
                    end
                        
                    
                    % Check that the two are the same to within 1%
                    diff_tolerance = 0.01;
                    rdel = reldiff(struct2array(ffit_check), struct2array(ffit));
                    if all(abs(rdel(fit_check_inds)) < diff_tolerance)
                        locs(a).fit_info = struct('ffit',ffit,'emgfit',emgfit,'param_stats',param_stats,'f0',f0,'history',history,'fitresults',fitresults);
                        break
                    elseif safety_count > max_fit_attempts
                        msg = sprintf('Could not fit %s in %d attempts', locs(a).ShortName, max_fit_attempts);
                        if fatal_if_cannot_fit
                            E.callError('fit_failure', msg);
                        else
                            fprintf('%s\n', msg);
                            locs(a).fit_info = misc_emissions_analysis.default_fit_structure;
                            break
                        end
                    else
                        fprintf('Attempt %d of %d: fit results differ by > %f%% (%s vs %s); retrying\n', safety_count, max_fit_attempts, diff_tolerance*100, struct2string(ffit), struct2string(ffit_check));
                        safety_count = safety_count + 1;
                    end
                end
                
                if ~isequal(locs(a).fit_info, misc_emissions_analysis.default_fit_structure)
                    % Add the emissions and lifetime. Use the 95% confidence
                    % intervals as the uncertainty. We need to restrict the
                    % winds to what should have been used to calculate the line
                    % densities.
                    [ ~, ~, ~, winds_strings ] = misc_emissions_analysis.extract_info_from_file_names(ldens_file);
                    % Assuming that "winds_strings" is of the form
                    % winds-lt# or winds-gt#, then
                    % winds_strings(end-2:end-1) will give the wind mode
                    % (less than or greater than) and converting
                    % winds_strings(3) to a number will give the speed.
                    wind_logical = misc_emissions_analysis.set_wind_conditions(locs(a), str2double(winds_strings(end)), winds_strings(end-2:end-1));
                    % We can use the WindUsedBool field if available to
                    % further constrain the winds to those from times when
                    % there were valid NO2 observations
                    if isfield(locs(a),'WindUsedBool')
                        wind_logical = wind_logical & locs(a).WindUsedBool;
                    else
                        warning('No "WindUsedBool" field detected; all winds that meet the speed criterion will be used in the lifetime/emissions calcuation');
                    end
                    emis_tau = misc_emissions_analysis.calculate_emission_lifetime(locs(a).no2_sectors, locs(a).fit_info, locs(a).WindSpeed(wind_logical));
                    
                    if add_nei
                        % Calculate the across-wind distance from the rotated
                        % latitude grid - since we rotate to the "x" (i.e.
                        % east-west) axis, across wind == latitudinally
                        across_wind_radius = abs((locs(a).no2_sectors.lat(1,1) - locs(a).no2_sectors.lat(end,1))/2);
                        
                        % Now get the WRF grid cells within that radius of the
                        % site and add up their NEI NO emissions.
                        xx = sqrt((nei_lon - locs(a).Longitude).^2 + (nei_lat - locs(a).Latitude).^2) < across_wind_radius;
                        
                        emis_tau.nei_emis = nansum(nei_avg_no(xx));
                    end
                    
                    locs(a).emis_tau = emis_tau;
                else
                    locs(a).emis_tau = misc_emissions_analysis.default_emis_tau_structure;
                end
            end
            
            dvec = line_densities.dvec;
            write_date = datestr(now);
            
            save(save_name, '-v7.3', 'locs', 'dvec', 'write_date');
        end
        
        function make_average_wrf_profiles(varargin)
            
            p = advInputParser;
            p.addParameter('time_period','');
            p.addParameter('relevant_hours', 15:23);
            p.addParameter('num_workers', nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            time_period = pout.time_period;
            hours_to_load = pout.relevant_hours;
            num_workers = pout.num_workers;
            active_pool = gcp('nocreate');
            if isnan(num_workers)
                if isempty(active_pool)
                    num_workers = 0;
                else
                    num_workers = active_pool.NumWorkers;
                end
            elseif isempty(active_pool)
                active_pool = parpool(num_workers);
            end
            
            if iscell(time_period)
                start_dates = time_period(1,:);
                end_dates = time_period(2,:);
            else
                [start_dates, end_dates] = misc_emissions_analysis.select_start_end_dates(time_period);
            end
            
            dvec = make_datevec(start_dates, end_dates);
            % This needs to be separate because we need a non-distributed
            % version for use outside the spmd block.
            distributed_dvec = distributed(dvec);
            profile_running_avgs = Composite(num_workers);
            vars_to_load = {'no','no2','hcho','pres'};
            n_vars = numel(vars_to_load);
            spmd(num_workers)
                local_dvec = getLocalPart(distributed_dvec);
                fprintf('Worker %d will do %s to %s\n', labindex, datestr(local_dvec(1)), datestr(local_dvec(end)));
                local_avgs = cell(1, n_vars);
                for i_var = 1:n_vars
                    local_avgs{i_var} = RunningAverage();
                end
                
                for i_date = 1:numel(local_dvec)
                    fprintf('  Loading files from %s\n', datestr(local_dvec(i_date)));
                    wrf_date = local_dvec(i_date);
                    wrf_dir = find_wrf_path('us','daily',wrf_date);
                    wrf_files = dir(fullfile(wrf_dir, sprintf('wrfout_*_%s*',datestr(wrf_date, 'yyyy-mm-dd'))));
                    wrf_hours = hour(date_from_wrf_filenames(wrf_files));
                    xx = ismember(wrf_hours, hours_to_load);
                    wrf_files(~xx) = [];
                    data = read_wrf_vars(wrf_dir, wrf_files, vars_to_load, 'squeeze', 0, 'as_struct');
                    if i_date == 1
                        xlon = ncread(fullfile(wrf_dir, wrf_files(1).name), 'XLONG');
                    end
                    
                    for i_var = 1:n_vars
                        tmp_var = wrf_day_weighted_average(xlon, 13.5, hours_to_load, data.(vars_to_load{i_var}));
                        local_avgs{i_var}.addData(tmp_var{1});
                    end
                end
                
                % Setting the value of a composite variable inside and SPMD
                % loop sets the value for this worker in the overall
                % composite.
                %
                % We want to pass the running averages back out because we
                % need to add up the running averages across all of the
                % workers.
                profile_running_avgs = local_avgs;
            end
            
            data = make_empty_struct_from_cell(vars_to_load);
            for i_var = 1:n_vars
                var_avg = RunningAverage();
                for i_worker = 1:numel(profile_running_avgs)
                    sub_avg = profile_running_avgs{i_worker};
                    % since adding data to a running average multiplies it
                    % by the weight, and getting the weighted average
                    % divides by the weight, to "append" weighted averages
                    % like this, we need to add the average, weighted by
                    % its weights, or we end up double-counting the weights
                    var_avg.addData(sub_avg{i_var}.getWeightedAverage(), sub_avg{i_var}.weights);
                end
                data.(vars_to_load{i_var}) = var_avg.getWeightedAverage();
            end
            
            % Add longitude and latitude
            wrf_file = find_wrf_path('us','daily',dvec(1),'fullpath');
            data.lon = ncread(wrf_file, 'XLONG');
            data.lat = ncread(wrf_file, 'XLAT');
            
            save(misc_emissions_analysis.wrf_avg_prof_file(start_dates, end_dates), '-struct', 'data', '-v7.3');
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % Plotting methods %
        %%%%%%%%%%%%%%%%%%%%
        
        function plot_site_summer_avg(varargin)
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('loc_to_plot', '');
            p.addParameter('plot_year',[]);
            p.addParameter('days_of_week', '');
            p.addParameter('monthly_or_daily', '');
            p.addParameter('plot_axis', gobjects(0));
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_to_plot = pout.loc_to_plot;
            plot_year = pout.plot_year;
            days_of_week = pout.days_of_week;
            monthly_or_daily = pout.monthly_or_daily;
            plot_axis = pout.plot_axis;
            
            locs = misc_emissions_analysis.read_locs_file();
            loc_names = {locs.ShortName};
            if isempty(loc_to_plot)
                loc_to_plot = ask_multichoice('Which location to plot?', loc_names, 'list', true);
            elseif ~any(strcmpi(loc_to_plot, loc_names))
                E.badinput('LOC_TO_PLOT must be one of the shortname in the trend locations file');
            end
            
            i_loc = strcmpi(loc_names, loc_to_plot);
            
            if isempty(plot_year)
                plot_year = ask_number('Enter the year (or years separated by a space) to do a summer average for', 'testfxn', @(x) all(x >= 2005 & x <= 2015), 'testmsg', 'Year(s) must be between 2005 and 2015');
            elseif ~isnumeric(plot_year)
                E.badinput('PLOT_YEAR must be numeric')
            end
            
            days_of_week = misc_emissions_analysis.choose_days_of_week(days_of_week);
            file_to_plot = misc_emissions_analysis.avg_file_name(plot_year, days_of_week);
            if ~exist(file_to_plot, 'file')
                E.badinput('No average file for year %d and days of week %s', plot_year, days_of_week);
            end   
            
            allowed_mod_strings = {'both', 'monthly', 'daily'};
            if isempty(monthly_or_daily)
                monthly_or_daily = 'both';
            elseif ~ismember(monthly_or_daily, allowed_mod_strings)
                E.badinput('MONTHLY_OR_DAILY must be one of: %s', allowed_mod_strings);
            end
            
            % From the average file, find the point near the given
            % location, then go out to 3x the radius given in the locs file
            avgs = load(file_to_plot);
            
            % Radius is in km, the lon and lats are in degrees. Assume ~110
            % km/deg.
            grid_del = abs(diff(avgs.monthly.lon(1,1:2)));
            radius_deg = locs(i_loc).Radius / 110;
            n_cells = ceil(radius_deg * 3 / grid_del);
            [yy, xx] = misc_emissions_analysis.find_indicies_in_box_around_point(locs(i_loc), avgs.monthly.lon, avgs.monthly.lat, n_cells);
            
            loc_longrid = avgs.monthly.lon(yy,xx);
            loc_latgrid = avgs.monthly.lat(yy,xx);
            if ismember(monthly_or_daily, {'both', 'monthly'})
                loc_no2grid = avgs.monthly.no2(yy,xx);
                figure;
                pcolor(loc_longrid, loc_latgrid, loc_no2grid);
                line(locs(i_loc).Longitude, locs(i_loc).Latitude, 'marker', 'p', 'color', 'k', 'linestyle', 'none');
                cb=colorbar;
                cb.Label.String = 'NO_2 VCD (molec. cm^{-2})';
                set(gca,'fontsize',16);
                shading flat;
                caxis(calc_plot_limits(loc_no2grid(:), 1e15, 'zero', 'max', [0 Inf]));
                title(sprintf('%s - Monthly profiles (%s, %s)', locs(i_loc).ShortName, sprintf_ranges(plot_year, 'value_sep', ', '), days_of_week));
            end
            % If daily profile avgs are available too, plot them as well
            if ~isscalar(avgs.daily.no2) && ismember(monthly_or_daily, {'both', 'daily'}) % if no data, the no2 grid will just be a scalar NaN
                loc_no2grid = avgs.daily.no2(yy,xx);
                if isempty(plot_axis)
                    figure; 
                    pcolor(loc_longrid, loc_latgrid, loc_no2grid);
                else
                    pcolor(plot_axis, loc_longrid, loc_latgrid, loc_no2grid);
                end
                line(locs(i_loc).Longitude, locs(i_loc).Latitude, 'marker', 'p', 'color', 'k', 'linestyle', 'none');
                cb=colorbar;
                cb.Label.String = 'NO_2 VCD (molec. cm^{-2})';
                set(gca,'fontsize',16);
                shading flat;
                caxis(calc_plot_limits(loc_no2grid(:), 1e15, 'zero', 'max', [0 Inf]));
                title(sprintf('%s - Daily profiles (%s %s)', locs(i_loc).ShortName, sprintf_ranges(plot_year, 'value_sep', ', '), days_of_week));
            end
        end
        
        function plot_both_sectors_interactive()
            % This will load the sites NO2 sectors file once then
            % continuously loop and let the user pick which location to
            % plot. It will then plot both the line densities and the
            % column densities by sector to help me decide which directions
            % to keep.
            avail_ld_files = dir(fullfile(misc_emissions_analysis.line_density_dir, '*.mat'));
                
            if numel(avail_ld_files) == 1
                ld_file = avail_ld_files(1).name;
            else
                ld_file = ask_multichoice('Select the sector line density file to use', {avail_ld_files.name}, 'list', true);
            end
            
            ld_file = fullfile(misc_emissions_analysis.line_density_dir, ld_file);
            locs_data = load(ld_file);
            locs_available = locs_data.locs;
            locs_dvec = locs_data.dvec;
            
            while true
                loc_inds = ask_multiselect('Choose the site(s) to plot', {locs_available.ShortName}, 'returnindex', true);
                locs_to_plot = locs_available(loc_inds);
                
                figs(1) = misc_emissions_analysis.plot_site_sectors_linedens(locs_to_plot, locs_dvec);
                figs(2) = misc_emissions_analysis.plot_sector_no2avg_with_boxes(locs_to_plot);
                
                tilefigs;
                
                if ~ask_yn('Plot another location?')
                    break
                else
                    close(figs)
                end
            end
        end
        
        function sectors_fig = plot_site_sectors_linedens(locs_to_plot, locs_dvec)
            E = JLLErrors;
            if nargin < 2
                avail_ld_files = dir(fullfile(misc_emissions_analysis.line_density_dir, '*.mat'));
                
                if numel(avail_ld_files) == 1
                    ld_file = avail_ld_files(1).name;
                else
                    ld_file = ask_multichoice('Select the sector line density file to use', {avail_ld_files.name}, 'list', true);
                end
                
                ld_file = fullfile(misc_emissions_analysis.line_density_dir, ld_file);
                locs_data = load(ld_file);
                locs_to_plot = locs_data.locs;
                locs_dvec = locs_data.dvec;
                
                loc_inds = ask_multiselect('Choose the site(s) to plot', {locs_to_plot.ShortName}, 'returnindex', true);
                locs_to_plot = locs_to_plot(loc_inds);
                
            elseif ~isstruct(locs_to_plot) || ~isfield(locs_to_plot, 'no2_sectors')
                E.badinput('LOCS_TO_PLOT must be a structure with field "no2_sectors"');
            end
            
            % Also load the summer average column density file so that we
            % can plot that as the center figure.
            locs_year = unique(year(locs_dvec));
            
            % Map the subplot index to the proper direction
            direction_names = {'NW','N','NE','W','vcds','E','SW','S','SE'};
            for a=1:numel(locs_to_plot)
                sectors_fig = figure;
                for b=1:9
                    ax=subplot(3,3,b);
                    if strcmpi(direction_names{b}, 'vcds')
                        misc_emissions_analysis.plot_site_summer_avg(locs_to_plot(a).ShortName, locs_year, 'daily', ax);
                        cb=colorbar;
                        cb.Label.String = 'NO_2 VCD (molec. cm^{-2})';
                        line(locs_to_plot(a).Longitude, locs_to_plot(a).Latitude, 'linestyle','none','marker','p','linewidth',2,'color','k');
                    else
                        plot(locs_to_plot(a).no2_sectors.x.(direction_names{b}), locs_to_plot(a).no2_sectors.linedens.(direction_names{b}));
                        xlabel('Dist. to site (km)');
                        ylabel('Line density (mol km^{-1})');
                        title(direction_names{b});
                    end
                end
            end
        end
        
        function sectors_fig = plot_sector_no2avg_with_boxes(locs_to_plot)
            if ~exist('locs_to_plot', 'var')
                avail_ld_files = dir(fullfile(misc_emissions_analysis.line_density_dir, '*.mat'));
                
                if numel(avail_ld_files) == 1
                    ld_file = avail_ld_files(1).name;
                else
                    ld_file = ask_multichoice('Select the sector line density file to use', {avail_ld_files.name}, 'list', true);
                end
                
                ld_file = fullfile(misc_emissions_analysis.line_density_dir, ld_file);
                locs_data = load(ld_file);
                locs_to_plot = locs_data.locs;
                
                loc_inds = ask_multiselect('Choose the site(s) to plot', {locs_to_plot.ShortName}, 'returnindex', true);
                locs_to_plot = locs_to_plot(loc_inds);
                
            elseif ~isstruct(locs_to_plot) || ~isfield(locs_to_plot, 'no2_sectors')
                E.badinput('LOCS_TO_PLOT must be a structure with field "no2_sectors"');
            end
            
            % Loop through the directions, plotting the average NO2 VCDs
            % with four different size boxes.
            direction_names = {'NW','N','NE','W','','E','SW','S','SE'};
            direction_angles = [135, 90, 45, 180, NaN, 0, -135, -90, -45];
            for a=1:numel(locs_to_plot)
                sectors_fig = figure;
                boxes_rel_x = [-0.5 1 1 -0.5 -0.5;...
                               -1 2 2 -1 -1;...
                               -2 4 4 -2 -2];
                boxes_rel_y = [-0.5 -0.5 0.5 0.5 -0.5;...
                               -1 -1 1 1 -1;...
                               -2 -2 2 2 -2];
                for b=1:9
                    subplot(3,3,b);
                    if isempty(direction_names{b})
                        title(locs_to_plot(a).Location);
                        axis off
                        continue
                    end
                    [lon, lat] = misc_emissions_analysis.rotate_lon_lat(locs_to_plot(a).no2_sectors.lon, locs_to_plot(a).no2_sectors.lat, locs_to_plot(a).Longitude, locs_to_plot(a).Latitude, direction_angles(b));
                    box_lon = boxes_rel_x + locs_to_plot(a).Longitude;
                    box_lat = boxes_rel_y + locs_to_plot(a).Latitude;
                    [box_lon, box_lat] = misc_emissions_analysis.rotate_lon_lat(box_lon, box_lat, locs_to_plot(a).Longitude, locs_to_plot(a).Latitude, direction_angles(b));
                    
                    pcolor(lon, lat, locs_to_plot(a).no2_sectors.no2_mean.(direction_names{b}));
                    shading flat; colorbar
                    for c=1:size(box_lon,1)
                        line(box_lon(c,:), box_lat(c,:), 'linewidth', 2, 'linestyle', '--', 'color', 'k');
                    end
                    title(direction_names{b});
                end
            end
        end
        
        function plot_sat_nei_emissions()
            % Plot a bar graph of OMI derived and NEI derived emissions.
            
            loc_inds = misc_emissions_analysis.get_loc_inds_interactive();
            
            if numel(loc_inds) > 1
                allowed_time_periods = {'beginning','end','both'};
                plot_time_period = ask_multichoice('Which time period(s) to plot?', allowed_time_periods, 'list', true);
                plot_beginning = any(strcmpi(plot_time_period, {'beginning','both'}));
                plot_end = any(strcmpi(plot_time_period, {'end','both'}));
                time_inds = [plot_beginning, plot_end];
            else
                time_inds = true(1,2);
            end
            
            
            
            [changes, loc_names] = misc_emissions_analysis.collect_changes('beginning','end','UMTWRFS','UMTWRFS','loc_inds',loc_inds);
            if numel(loc_inds) == 1
                % If only plotting one location, then we can split up the
                % bars more nicely than if plotting multiple locations
                plot_data = cat(2, changes.emis', changes.nei_emis');
                legend_str = {'Top-down (BEHR)','Bottom-up (NEI)'};
                legend_inds = true(size(legend_str));
                xticklabels =  {'2005,07','2012-13'};
                axis_opts = {'xticklabels', xticklabels(time_inds), 'fontsize',14};
            else
                plot_data = cat(2, changes.emis(:,time_inds), changes.nei_emis(:,time_inds));
                legend_inds = [time_inds(1), time_inds(2), time_inds(1), time_inds(2)]; % needed for the legend string
                legend_str = {'05-07 Top-down (BEHR)','12-13 Top-down (BEHR)','05-07 Bottom-up (NEI)','12-13 Bottom-up (NEI)'};
                axis_opts = {'xticklabels', loc_names', 'xticklabelrotation', 30, 'fontsize', 14};
            end
            
            figure;
            bar(plot_data);
            set(gca, axis_opts{:});
            legend(legend_str{legend_inds});
            ylabel('NO Emissions (Mg h^{-1})');
            
        end
        
        function plot_fits_interactive(varargin)
            p = inputParser;
            p.addParameter('loc_indicies', [])
            p.addParameter('include_2yr', nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            loc_inds = pout.loc_indicies;
            include_2yr = opt_ask_yn('Include 2yr periods?', pout.include_2yr, '"include_2yr"');
            
            % Load the weekday, weekend, and all day fits files. Cut out
            % the memory intensive parts of the line density sub structure
            % to save memory.
            [use_wrf, loc_inds_wrf] = misc_emissions_analysis.ask_to_use_wrf();
            if isempty(loc_inds)
                loc_inds = loc_inds_wrf;
            end
            
            [periods.beg.start_date, periods.beg.end_date] = misc_emissions_analysis.select_start_end_dates('beginning');
            periods.beg.loc_inds = loc_inds;
            [periods.end.start_date, periods.end.end_date] = misc_emissions_analysis.select_start_end_dates('end');
            periods.end.loc_inds = loc_inds;
            if include_2yr
                [periods.beg2yr.start_date, periods.beg2yr.end_date] = misc_emissions_analysis.select_start_end_dates('beg_2yr');
                periods.beg2yr.loc_inds = 1:70;  % The 2-year files only have 70 locations
                [periods.end2yr.start_date, periods.end2yr.end_date] = misc_emissions_analysis.select_start_end_dates('end_2yr');
                periods.end2yr.loc_inds = 1:70;
            end
            
            fit_type = misc_emissions_analysis.get_fit_type_interactive();
            
            
            
            periods_fns = fieldnames(periods);
            for i_period = 1:numel(periods_fns)
                per_fn = periods_fns{i_period};
                days_of_week = {'UMTWRFS', 'TWRF', 'US'};
                loc_prototype = struct('Location', '', 'x', [], 'linedens', [], 'emgfit', [], 'r2', []);
                periods.(per_fn).locs = repmat(loc_prototype, numel(loc_inds), numel(days_of_week));
                
                for i_dow = 1:numel(days_of_week)
                    the_fits = load(misc_emissions_analysis.fits_file_name(periods.(per_fn).start_date, periods.(per_fn).end_date, use_wrf, periods.(per_fn).loc_inds, days_of_week{i_dow}, fit_type));
                    for i_loc = 1:numel(the_fits.locs)
                        periods.(per_fn).locs(i_loc, i_dow) = copy_structure_fields(the_fits.locs(i_loc), loc_prototype, 'substructs');
                    end
                end
            end
            % Now we can actually do the plotting. Make a list of the
            % available locations, then ask which one to plot until we quit
            locs_list = {periods.(periods_fns{1}).locs(:,1).Location};
            colors = {[0.5 0.5 0.5], 'c', [1 0.5 0];...
                      'k'          , 'b', 'r'};
            while true
                loc_ind = ask_multichoice('Plot which location?', locs_list, 'list', true, 'index', true);
                figs = gobjects(size(periods_fns));
                for i_period = 1:numel(periods_fns)
                    figs(i_period) = figure;
                    this_period = periods.(periods_fns{i_period});
                    l = gobjects(size(this_period.locs,2),1);
                    for i_wkday = 1:size(this_period.locs,2)
                        l(i_wkday) = line(this_period.locs(loc_ind, i_wkday).x, this_period.locs(loc_ind, i_wkday).linedens, 'color', colors{1, i_wkday}, 'marker', 'o', 'linestyle', 'none');
                        if ~isempty(this_period.locs(loc_ind, i_wkday).emgfit)
                            line(this_period.locs(loc_ind, i_wkday).x, this_period.locs(loc_ind, i_wkday).emgfit, 'color', colors{2, i_wkday}, 'linestyle', '--');
                        end
                    end
                    r2_array = {this_period.locs(loc_ind, :).r2};
                    xx = iscellcontents(r2_array, @isempty);
                    r2_array(xx) = {nan};
                    legend(l, sprintfmulti('%s (R^2 = %.2f)', days_of_week, r2_array));
                    title(sprintf('%s, %s', periods_fns{i_period}, this_period.locs(loc_ind, 1).Location));
                end
                input('Press ENTER to continue','s');
                
                close(figs);
            end
        end
        
        function figs = plot_plume_differences(plot_type, varargin)
            E = JLLErrors;
            
            p = advInputParser;
            p.addParameter('time_period_1', '');
            p.addParameter('time_period_2', '');
            p.addParameter('days_of_week_1', '');
            p.addParameter('days_of_week_2', '');
            p.addParameter('use_wrf_1', nan);
            p.addParameter('use_wrf_2', nan);
            p.addParameter('loc_indicies', nan);
            p.addParameter('separate_figs', false);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            separate_figs = pout.separate_figs;
            
            loc_inds = misc_emissions_analysis.ask_for_loc_inds(pout.loc_indicies);
            
            [start_1, end_1, ~, tp_1_legend] = misc_emissions_analysis.select_start_end_dates(pout.time_period_1, 'Select the time period for the first plume.');
            [days_of_week_1, dow_1_legend] = misc_emissions_analysis.select_days_of_week(pout.days_of_week_1, 'Select the days of week for the first plume.');
            use_wrf_1 = opt_ask_yn('Use WRF for the first plume? (BEHR if not)', pout.use_wrf_1);
            tp_1_legend = sprintf('%s %s', tp_1_legend, dow_1_legend);
            
            [start_2, end_2, ~, tp_2_legend] = misc_emissions_analysis.select_start_end_dates(pout.time_period_2, 'Select the time period for the second plume.');
            [days_of_week_2, dow_2_legend] = misc_emissions_analysis.select_days_of_week(pout.days_of_week_2, 'Select the days of week for the second plume.');
            use_wrf_2 = opt_ask_yn('Use WRF for the second plume? (BEHR if not)', pout.use_wrf_2);
            tp_2_legend = sprintf('%s %s', tp_2_legend, dow_2_legend);
            
            % I made symlinks to the 2-year periods that use locs1-71 in
            % the name to make this easier.
            file1 = misc_emissions_analysis.fits_file_name(start_1, end_1, use_wrf_1, 1:71, days_of_week_1, 'lu');
            file2 = misc_emissions_analysis.fits_file_name(start_2, end_2, use_wrf_2, 1:71, days_of_week_2, 'lu');
            
            base = load(file1);
            new = load(file2);
            
            base.locs = misc_emissions_analysis.cutdown_locs_by_index(base.locs, loc_inds);
            new.locs = misc_emissions_analysis.cutdown_locs_by_index(new.locs, loc_inds);
            
            switch lower(plot_type)
                case '2d'
                    figs = misc_emissions_analysis.plot_plume_difference_2d(base, new, separate_figs);
                case '1d'
                    figs = misc_emissions_analysis.plot_plume_difference_1d(base, tp_1_legend, new, tp_2_legend, separate_figs);
                otherwise
                    E.notimplemented('No plotting action implemented for "%s"', plot_type)
            end
            
        end
        
        function figs = plot_plume_difference_1d(base, base_legstr, new, new_legstr, separate_figs)
            E = JLLErrors;
            n_locs = numel(base.locs);
            if n_locs ~= numel(new.locs)
                E.badinput('BASE and NEW must have the same number of locations')
            end
            
            if separate_figs
                figs = gobjects(1, n_locs);
            else
                figs = figure;
            end
            
            marker = 'o';
            marker_size = 6;
            marker_line_width = 1.5;
            fit_line_width = 2;
            
            for i_loc = 1:n_locs
                next_plot();
                
                base_x = base.locs(i_loc).no2_sectors.x;
                base_ld = range_squeeze(base.locs(i_loc).no2_sectors.linedens, [0 1]);
                base_fit = range_squeeze(base.locs(i_loc).fit_info.emgfit, [0 1]);
                new_x = base.locs(i_loc).no2_sectors.x;
                new_ld = range_squeeze(new.locs(i_loc).no2_sectors.linedens, [0 1]);
                new_fit = range_squeeze(new.locs(i_loc).fit_info.emgfit, [0 1]);
                
                line(base_x, base_ld, 'color', [0 0.5 0], 'marker', marker, 'linewidth', marker_line_width, 'markersize', marker_size, 'linestyle', 'none');
                line(base_x, base_fit, 'color', 'b', 'marker', 'none', 'linewidth', fit_line_width);
                line(new_x, new_ld, 'color', [1 0.5 0], 'marker', marker, 'linewidth', marker_line_width, 'markersize', marker_size, 'linestyle', 'none');
                line(new_x, new_fit, 'color', 'r', 'marker', 'none', 'linewidth', fit_line_width);
                
                legend(cprintf({'Line density (%s)', 'Fit (%s)', 'Line density (%s)', 'Fit (%s)'}, {base_legstr, base_legstr, new_legstr, new_legstr}));
                set(gca,'fontsize', 12);
                ylabel('Normalized NO_2');
                xlabel('Distance from city (km)');
            end
            
            
            function ax = next_plot()
                if separate_figs
                    figs(i_loc) = figure;
                    ax = gca;
                else
                    ax = subplot(1, n_locs, i_loc);
                end
            end
        end
        
        function figs = plot_plume_difference_2d(base, new, separate_figs)
            % Not intended to be called directly, but you can if you want.
            % Plots the differences in the 2D plumes between the base and
            % new times for each location in them. BASE and NEW are the
            % structure resulting from loading one of the fits files,
            % SEPARATE_FIGS is a boolean if you want each plot as a
            % separate figure or all on one plot.
            
            E = JLLErrors;
            n_locs = numel(base.locs);
            if n_locs ~= numel(new.locs)
                E.badinput('BASE and NEW must have the same number of locations')
            end
            
            if separate_figs
                figs = gobjects(3, n_locs);
            else
                figs = figure;
            end
            
            abs_colormap = jet;
            diff_colormap = blue_red_cmap;
            
            for i_loc = 1:n_locs
                % Calculate the best limits for the absolute value plots
                abs_clims = calc_plot_limits(veccat(base.locs(i_loc).no2_sectors.no2_mean(:), new.locs(i_loc).no2_sectors.no2_mean(:)), 1e15, 'zero');
                
                curr_ax = next_plot(1);
                plot_plume(curr_ax, base.locs(i_loc));
                
                curr_ax = next_plot(2);
                plot_plume(curr_ax, new.locs(i_loc));
                
            end
            
            if ~separate_figs
                label_subfigs(figs, 'xshift', 0.2);
            end
            
            function ax = next_plot(y)
                if separate_figs
                    figs(y,i_loc) = figure;
                    ax = gca;
                else
                    idx = (y-1)*n_locs + i_loc;
                    ax = subplot(2, n_locs, idx);
                end
            end
            
            function plot_plume(ax, loc, second_loc)
                if nargin > 2
                    if ~isequal(second_loc.no2_sectors.lon, loc.no2_sectors.lon) || ~isequal(second_loc.no2_sectors.lat, loc.no2_sectors.lat)
                        E.badinput('loc and second_loc must have the same lat/lon')
                    end
                    
                    plot_quantity = second_loc.no2_sectors.no2_mean - loc.no2_sectors.no2_mean;
                    clims = calc_plot_limits(plot_quantity, 1e15, 'diff');
                    cmap = diff_colormap;
                    cb_label = '\Delta VCD (molec cm^{-2})';
                else
                    plot_quantity = loc.no2_sectors.no2_mean;
                    clims = abs_clims;
                    cmap = abs_colormap;
                    cb_label = 'VCD (molec cm^{-2})';
                end
                pcolor(ax, loc.no2_sectors.lon, loc.no2_sectors.lat, plot_quantity);
                line(loc.Longitude, loc.Latitude, 'linestyle', 'none', 'marker', 'p', 'color', 'k', 'markerfacecolor','k','markersize',16);
                shading flat;
                colormap(ax, cmap);
                caxis(clims);
                cb = colorbar;
                cb.Label.String = cb_label;
                
                % Also find, from the line densities where it drops to 1/e
                % of the max value. Only do this where plotting a single
                % plume, not a difference.
                
                x = loc.no2_sectors.x;
                linedens = loc.no2_sectors.linedens;
                linedens_boxcar_mean = movmean(linedens, 5);
                
                [max_ld, i_max] = max(linedens);
                background_ld = loc.fit_info.ffit.B;
                % Find the distance where the boxcar average first brackets
                % the 1/e value. We'll interpolate then to get the exact
                % distance.
                i_x = 1:numel(linedens);
                % We want the line density where the enhancement has
                % dropped to 1/e of its original value, so we caculate 1/e
                % of that and add it back into the background.
                e_folding_ld = exp(-1) * (max_ld - background_ld) + background_ld;
                e_idx = find(linedens_boxcar_mean < e_folding_ld & i_x > i_max, 1, 'first');
                x_e_fold = interp1(linedens_boxcar_mean(e_idx-1:e_idx), x(e_idx-1:e_idx), e_folding_ld);
                
                % Now we need to draw a vertical line on the plot at the
                % e-folding distance downwind (maybe should draw a line
                % between the maximum and the e-folding distance?)
                % First step is to convert the distances back to
                % longitudes. 
                
                % m_fdist expects distances in meters, hence multiply x by
                % 1000. m_fdist also returns longitude in degrees east (so
                % [0, 360] instead of [-180, 180], so we need to convert
                % that back.
                lon_max = m_fdist(loc.Longitude, loc.Latitude, 90, x(i_max)*1000) - 360;
                lon_efold = m_fdist(loc.Longitude, loc.Latitude, 90, x_e_fold*1000) - 360;
                
                % We'll put the line halfway towards the bottom of the plot
                bottom_lat = min(ax.YLim);
                line_lat = mean([bottom_lat, loc.Latitude]);
                
                line_grp = hggroup(ax);
                line_args = {'color', 'w', 'linewidth', 4, 'parent', line_grp};
                inner_line_args = {'color', 'k', 'linewidth', 2, 'parent', line_grp};
                axis_line_args = {'color', 'w', 'linewidth', 2, 'linestyle', '--', 'parent', line_grp};
                
                % Also draw dashed lines down to the x-axis so people can
                % compare the e-folding absolute distances more easily.
                % Draw these first so they end up on the bottom
                line([lon_max, lon_max], [line_lat, bottom_lat], axis_line_args{:});
                line([lon_efold, lon_efold], [line_lat, bottom_lat], axis_line_args{:});
                
                line([lon_max, lon_efold], [line_lat, line_lat], line_args{:});
                line([lon_max, lon_efold], [line_lat, line_lat], inner_line_args{:});
                % Make tips on the end
                tip_radius = diff(ax.YLim)/40;
                line([lon_max, lon_max], [line_lat - tip_radius, line_lat + tip_radius], line_args{:});
                line([lon_max, lon_max], [line_lat - tip_radius, line_lat + tip_radius], inner_line_args{:});
                line([lon_efold, lon_efold], [line_lat - tip_radius, line_lat + tip_radius], line_args{:});
                line([lon_efold, lon_efold], [line_lat - tip_radius, line_lat + tip_radius], inner_line_args{:});
                
                set(ax, 'fontsize', 14);
            end
        end
        
        function plot_emis_tau_changes(varargin)
            
            [first_time_period, second_time_period, first_weekdays, second_weekdays, loc_types, series_labels] = misc_emissions_analysis.get_change_file_selection_input(varargin{:});
            
            [changes, loc_names] = misc_emissions_analysis.collect_changes(first_time_period, second_time_period, first_weekdays, second_weekdays, 'loc_types', loc_types);
            
            plot_grouped_changes({changes.emis, changes.nei_emis}, 'group_labels', {'BEHR','NEI'}, 'series_labels', series_labels, 'inter_space', 2, 'tick_labels', loc_names);
            set(gca,'XTickLabelRotation',30);
            plot_grouped_changes({changes.tau}, 'group_labels', {'Lifetime'}, 'series_labels', series_labels, 'inter_space', 2, 'tick_labels', loc_names);
            set(gca,'XTickLabelRotation',30);
        end
        
        function plot_emis_change_map(varargin)
            default_differences = struct('emis_type', {'emis','nei_emis';'emis','nei_emis'},...
                                         'time_period', {'beginning', 'beginning'; 'end', 'end'},...
                                         'days_of_week', {'UMTWRFS', 'UMTWRFS'; 'UMTWRFS', 'UMTWRFS'});
            
            p = inputParser;
            p.addParameter('differences',default_differences);
            p.addParameter('loc_types','');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_types = pout.loc_types;
            differences = pout.differences;
            
            allowed_loc_types = {'Cities','PowerPlants'};
            if ~ischar(loc_types) && ~iscellstr(loc_types)
                E.badinput('"loc_types" must be a character array or cell array of character arrays')
            elseif ~isempty(loc_types)
                if ~all(ismember(loc_types, allowed_loc_types))
                    E.badinput('"loc_types" must be one or more of the following: %s', strjoin(allowed_loc_types, ', '));
                end
            else
                loc_types = ask_multiselect('Select one or more site types to plot:', allowed_loc_types);
            end
            
            n_differences = size(differences, 1);
            plot_loc_inds = misc_emissions_analysis.loc_types_to_inds(loc_types{:});
            % Add an extra point to the corners so that the last point
            % we actually plot does not overlap the first one, but all
            % points are evenly spaced around the circle.
            plot_corners = fliplr(linspace(-180,180,n_differences+1));
            
            map_fig=figure;
            map_ax = gca;
            state_outlines('k','not','ak','hi');
            max_diff = 0;
            for i_change = 1:n_differences
                [this_change, loc_names, loc_coords] = misc_emissions_analysis.collect_changes(differences(i_change,1).time_period, differences(i_change,2).time_period, differences(i_change,1).days_of_week, differences(i_change,2).days_of_week, 'loc_inds', plot_loc_inds);
                % for now, assume no uncertainty in the NEI emissions.
                this_change.nei_emis_sd = zeros(size(this_change.emis_sd));
                
                emis_field_1 = differences(i_change,1).emis_type;
                sd_field_1 = sprintf('%s_sd', emis_field_1);
                emis_field_2 = differences(i_change,2).emis_type;
                sd_field_2 = sprintf('%s_sd', emis_field_2);
                
                emis_values = [this_change.(emis_field_1)(:,1), this_change.(emis_field_2)(:,2)];
                sd_values = [this_change.(sd_field_1)(:,1), this_change.(sd_field_2)(:,2)];

                max_diff = ceil(max(max_diff, max(abs(diff(emis_values, [], 2)))));
                
                this_series = misc_emissions_analysis.create_map_series(emis_values, sd_values, this_change.n_dofs, this_change.r2, loc_coords, plot_corners(i_change), loc_names);
                if plot_corners(i_change) == min(abs(plot_corners))
                    % Only include names on the right most point.
                    misc_emissions_analysis.plot_map_series(map_ax, this_series, 'include_names');
                else
                    misc_emissions_analysis.plot_map_series(map_ax, this_series);
                end
                
                % Make box plots as well
                figure;
                boxplot(diff(emis_values,[],2));
                set(gca,'fontsize',16,'xticklabels',{''});
                title_emis_types = struct('emis', 'BEHR', 'nei_emis', 'NEI');
                title(sprintf('%s (%s, %s) - %s (%s, %s)', title_emis_types.(emis_field_2), differences(i_change,2).time_period, differences(i_change,2).days_of_week, title_emis_types.(emis_field_1), differences(i_change,1).time_period, differences(i_change,1).days_of_week));
                ylabel('\Delta E_{NO2} (Mg NO_x h^{-1})');
            end
            
            figure(map_fig);
            cb = colorbar;
            cb.Label.String = '\Delta E_{NO2} (Mg NO_x h^{-1}, NEI - BEHR)';
            colormap(blue_red_only_cmap);
            caxis([-max_diff max_diff]);
            set(gca,'fontsize',16);
        end
        
        function plot_lifetime_change_map(varargin)
            % Makes a scatter plot of changes in lifetime, both 2012-2013
            % minus 2005-2007 and weekend minus weekday, plotted on a map.
            % Each location will have three symbols in a triangle, the top
            % symbol will represent the decadal change, the bottom left one
            % the 2005-2007 weekend/weekday change, and the bottom right
            % right the 2012-2013 weekend/weekday change. Different symbols
            % will be used to represent the significance of each change;
            % circles for changes significant using 95% CIs, asterisks for
            % changes not significant using 95% CIs, and X's for changes in
            % which one of the R2 values is less than some critical value.
            
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('loc_types','');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_types = pout.loc_types;
            
            allowed_loc_types = {'Cities','PowerPlants'};
            if ~ischar(loc_types) && ~iscellstr(loc_types)
                E.badinput('"loc_types" must be a character array or cell array of character arrays')
            elseif ~isempty(loc_types)
                if ~all(ismember(loc_types, allowed_loc_types))
                    E.badinput('"loc_types" must be one or more of the following: %s', strjoin(allowed_loc_types, ', '));
                end
            else
                loc_types = ask_multiselect('Select one or more site types to plot:', allowed_loc_types);
            end
            
            error('Change in behavior: Need to convert loc_types to location indicies');
            
            weekdays = 'TWRF';
            weekends = 'US';
            alldays = 'UMTWRFS';
            
            [decadal_changes, loc_names, loc_coords] = misc_emissions_analysis.collect_changes('beginning', 'end', alldays, alldays, 'loc_types', loc_types);
            beginning_changes = misc_emissions_analysis.collect_changes('beginning', 'beginning', weekdays, weekends, 'loc_types', loc_types);
            end_changes = misc_emissions_analysis.collect_changes('end', 'end', weekdays, weekends, 'loc_types', loc_types);
            
            decadal_series = misc_emissions_analysis.create_map_series(decadal_changes.tau, decadal_changes.tau_sd, decadal_changes.n_dofs, decadal_changes.r2, loc_coords, 'top', loc_names);
            beginning_series = misc_emissions_analysis.create_map_series(beginning_changes.tau, beginning_changes.tau_sd, beginning_changes.n_dofs, beginning_changes.r2, loc_coords, 'bottom-left', loc_names);
            end_series = misc_emissions_analysis.create_map_series(end_changes.tau, end_changes.tau_sd, end_changes.n_dofs, end_changes.r2, loc_coords, 'bottom-right', loc_names);
            
            figure;
            state_outlines('k','not','ak','hi');
            misc_emissions_analysis.plot_map_series(gca, decadal_series, 'include_names');
            misc_emissions_analysis.plot_map_series(gca, beginning_series);
            misc_emissions_analysis.plot_map_series(gca, end_series);
            cb = colorbar;
            cb.Label.String = '\Delta \tau_{NO2} (hours)';
            colormap(blue_red_only_cmap);
            caxis([-6 6]);
            set(gca,'fontsize',16);
        end
        
        function figs = plot_lifetime_vs_mass(varargin)
            % Plot lifetimes vs. some measure of NOx mass for each location
            % separately. Parameters:
            %   'loc_inds' - numeric indicies of which locations to
            %   include.
            %
            %   'mass_value' - character array, either 'a' or 'vcds'.
            %   'a' uses the fitting parameter a, 'vcds' uses the
            %   average summer columns within the box width of the site.
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('loc_inds', nan);
            p.addParameter('mass_value', '');
            p.addParameter('fit_type', '');
            p.addParameter('sat_or_model', {});
            p.addParameter('single_plot', nan);
            p.addParameter('single_plot_mode', '');
            p.addParameter('include_2years', nan);
            p.addParameter('days_of_week', '');
            p.addParameter('connect_wkend', nan);
            p.addParameter('bad_fit_display','');
            p.addParameter('legend', '');
            p.addParameter('title', true); % if plotting interactively, we want the title. But provide the option to turn in off if not plotting interactively
            
            p.parse(varargin{:});
            pout = p.Results;
            
            include_title = pout.title;
            
            loc_inds = pout.loc_inds;
            if isnan(loc_inds)
                [loc_inds, file_loc_inds] = misc_emissions_analysis.get_loc_inds_interactive();
            else
                file_loc_inds = 1:71;
            end
            
            mass_value = pout.mass_value;
            allowed_mass_vals = {'a','a_plus_B','vcds'};
            if isempty(mass_value)
                mass_value = ask_multichoice('Which quantity to use for mass of NOx?', allowed_mass_vals, 'list', true);
            elseif ~ismember(mass_value, allowed_mass_vals)
                E.badinput('MASS_VALUE must be one of: %s', strjoin(allowed_mass_vals, ', '));
            end
            
            include_2years = opt_ask_yn('Include the points from the 2 year runs?', pout.include_2years, '"include_2years"');
            
            sat_or_model = opt_ask_multiselect('Which data source to use?', {'BEHR', 'WRF'}, pout.sat_or_model, '"sat_or_model"');
            
            single_plot_bool = opt_ask_yn('Plot all locations on a single plot?', pout.single_plot, '"single_plot"');
            if single_plot_bool
                single_plot_mode = opt_ask_multichoice('What should color represent in the plot?', {'Year','HCHO VCD'}, pout.single_plot_mode, '"single_plot_mode"', 'list', true);
            else
                single_plot_mode = '';
            end
            
            allowed_dows = {'UMTWRFS','TWRF','US'};
            days_of_week = opt_ask_multiselect('Choose which day-of-week subsets to include', [allowed_dows, 'all'], pout.days_of_week, '"days_of_week"');
            if strcmpi(days_of_week, 'all')
                days_of_week = allowed_dows;
            elseif ischar(days_of_week)
                days_of_week = {days_of_week};
            end
            
            do_connect_wkday_wkend = false;
            include_decade = ismember('UMTWRFS', days_of_week);
            include_wkday_wkend = all(ismember({'TWRF','US'}, days_of_week));
            include_behr = ismember('BEHR', sat_or_model);
            include_wrf = ismember('WRF', sat_or_model);
            if ~single_plot_bool
                if include_wkday_wkend
                    do_connect_wkday_wkend = opt_ask_yn('Connect weekday/weekend points?', pout.connect_wkend, '"connect_wkend"');
                elseif xor(ismember('TWRF', days_of_week), ismember('US', days_of_week))
                    E.notimplemented('For non-single plot mode, having one but not both of "TWRF" and "US" are not supported')
                end
                
                allowed_legend_figs = {'all', 'none', 'first', 'last'};
            else
                allowed_legend_figs = {'all', 'none'};
            end
            where_to_put_legend = opt_ask_multichoice('Where to include a legend?', allowed_legend_figs, pout.legend, '"legend"', 'list', true);
            
            bad_fit_display = opt_ask_multichoice('Plot bad fits?', {'no','grey'}, pout.bad_fit_display, '"bad_fit_display"', 'list', true');
            
            fit_type = misc_emissions_analysis.get_fit_type_interactive(pout.fit_type);
            
            locs = misc_emissions_analysis.read_locs_file();
            
            vcds_bool = strcmpi(mass_value, 'vcds');
            [~,~,~,wrf_file_inds] = misc_emissions_analysis.ask_to_use_wrf(true);
            
            % End input section %
            
            % These help keep the style consistent %
            marker_size = 10;
            marker_linewidth = 1.5;
            connector_linewidth = 2;
            dow_markers = misc_emissions_analysis.dow_markers;
            product_series = struct('behr', struct('used', false, 'name', 'BEHR', 'style', struct('marker', 'o', 'color', 'k', 'markerfacecolor', 'k', 'linestyle', 'none', 'linewidth', 2)),...
                'wrf', struct('used', false, 'name', 'WRF', 'style', struct('marker', 'o', 'color', 'k', 'linestyle','none','linewidth',2)));
            
            time_period_colors = misc_emissions_analysis.time_period_colors;
            % Now actually load the data
            
            all_changes = struct([]);
            if include_behr
                if include_decade
                    if include_2years
                        load_change_group('beg_2yr', 'end_2yr', 'UMTWRFS', 'UMTWRFS', 1:70);
                    end
                    
                    load_change_group('beginning', 'end', 'UMTWRFS', 'UMTWRFS');
                end
                if include_wkday_wkend
                    if include_2years
                        load_change_group('beg_2yr', 'beg_2yr', 'TWRF', 'US', 1:70);
                        load_change_group('end_2yr', 'end_2yr', 'TWRF', 'US', 1:70);
                    end
                    load_change_group('beginning', 'beginning', 'TWRF', 'US');
                    load_change_group('end', 'end', 'TWRF', 'US');
                end
            end
            
            if include_wrf
                % WRF is only processed for all days because it cannot have
                % a weekend effect, since the NEI emissions are just a 24
                % hour cycle
                load_change_group('beg_2yr', 'end_2yr', 'UMTWRFS', 'UMTWRFS', wrf_file_inds, true);
                load_change_group('beginning', 'end', 'UMTWRFS', 'UMTWRFS', wrf_file_inds, true);
            end
            
            if do_connect_wkday_wkend
                conn_fmt_fxn = @make_connector_fmt;
            else
                conn_fmt_fxn = @(x, y) struct('linestyle', 'none');
            end
            
            if vcds_bool
                x_label_str = 'Avg. NO_2 VCD (molec. cm^2)';
            else
                x_label_str = 'a (mol NO_2)';
            end
            
            if ~single_plot_bool
                figs = gobjects(numel(loc_inds),1);
                for i_loc = 1:numel(loc_inds)
                    figs(i_loc) = figure;
                    ax = gca;
                    for i_change = 1:numel(all_changes)
                        this_group_style = make_group_fmt(all_changes(i_change).style, all_changes(i_change).is_fit_good(i_loc, :));
                        this_connector_style = conn_fmt_fxn(all_changes(i_change).style, all_changes(i_change).is_significant(i_loc), all_changes(i_change).is_fit_good(i_loc,:));
                        plot_changes(all_changes(i_change).(mass_value)(i_loc,:), all_changes(i_change).tau(i_loc,:),...
                            'group_fmts', this_group_style, 'connector_fmt', this_connector_style, 'parent', ax);
                    end
                    
                    
                    if strcmpi(where_to_put_legend, 'all') || (strcmpi(where_to_put_legend, 'first') && i_loc == 1) || (strcmpi(where_to_put_legend, 'last') && i_loc == numel(loc_inds))
                        make_legend(ax);
                    end
                    if include_title
                        title(locs(loc_inds(i_loc)).Location);
                    end
                    
                    set(ax,'fontsize',16);
                    xlabel(x_label_str);
                    ylabel('\tau (hours)');
                end
            else
                l = gobjects(0);
                legend_cell = {};
                figure;
                if strcmpi(single_plot_mode, 'Year')
                    for i_change = 1:numel(all_changes)
                        line(all_changes(i_change).(mass_value)(:,1), all_changes(i_change).tau(:,1), all_changes(i_change).style(1));
                        line(all_changes(i_change).(mass_value)(:,2), all_changes(i_change).tau(:,2), all_changes(i_change).style(2));
                    end
                    legend_flags = {};
                elseif strcmpi(single_plot_mode, 'HCHO VCD')
                    is_wrf = [all_changes.is_wrf];
                    behr_no2_vcds = [all_changes(~is_wrf).(mass_value)];
                    behr_hcho_vcds = [all_changes(~is_wrf).hcho_vcds];
                    behr_tau = [all_changes(~is_wrf).tau];
                    scatter(behr_no2_vcds(:), behr_tau(:), 60, behr_hcho_vcds(:), 'filled');
                    
                    hold on
                    wrf_no2_vcds = [all_changes(is_wrf).(mass_value)];
                    wrf_hcho_vcds = [all_changes(is_wrf).hcho_vcds];
                    wrf_tau = [all_changes(is_wrf).tau];
                    scatter(wrf_no2_vcds(:), wrf_tau(:), 60, wrf_hcho_vcds(:));
                    
                    cb = colorbar;
                    cb.Label.String = 'HCHO VCD (molec. cm^{-2})';
                    colormap jet
                    legend_flags = {'no_dow', 'no_time_period'};
                end
                
                
                
                if strcmpi(where_to_put_legend, 'all')
                    make_legend(gca, legend_flags{:});
                end
                set(gca,'xscale','log','fontsize',14);
                xlabel(x_label_str);
                ylabel('\tau (hours)');
            end
            
            %                                                     %
            % Internal helper functions for plot_lifetime_vs_mass %
            %                                                     %
            
            function load_change_group(time_period_1, time_period_2, days_of_week_1, days_of_week_2, varargin)
                
                if numel(varargin) >= 1
                    load_loc_inds = varargin{1};
                else
                    load_loc_inds = 1:71;
                end
                
                if numel(varargin) >= 2
                    load_wrf = varargin{2};
                else
                    load_wrf = false;
                end
                
                if load_wrf
                    marker_fills = 'none';
                else
                    marker_fills = {time_period_colors.(time_period_1).color, time_period_colors.(time_period_2).color};
                end
                
                changes = misc_emissions_analysis.collect_changes(time_period_1, time_period_2, days_of_week_1, days_of_week_2, 'loc_inds', loc_inds, 'file_loc_inds', load_loc_inds, 'use_wrf', load_wrf, 'include_vcds', vcds_bool, 'fit_type', fit_type);
                changes.is_significant = misc_emissions_analysis.is_change_significant(changes.tau, changes.tau_sd, changes.n_dofs);
                changes.style = struct('marker', {dow_markers.(days_of_week_1).marker, dow_markers.(days_of_week_2).marker},...
                    'linestyle', 'none', 'color', {time_period_colors.(time_period_1).color, time_period_colors.(time_period_2).color},...
                    'markersize',marker_size,'linewidth',marker_linewidth,'markerfacecolor',marker_fills);
                changes.is_wrf = load_wrf;
                dow_markers.(days_of_week_1).used = true;
                dow_markers.(days_of_week_2).used = true;
                if load_wrf
                    product_series.wrf.used = true;
                else
                    product_series.behr.used = true;
                end
                time_period_colors.(time_period_1).used = true;
                time_period_colors.(time_period_2).used = true;
                
                
                if isempty(all_changes)
                    all_changes = changes;
                else
                    all_changes(end+1) = changes;
                end
            end
            
            function make_legend(parent, varargin)
                subp = advInputParser;
                subp.addFlag('no_dow');
                subp.addFlag('no_time_period');
                subp.addFlag('no_product');
                
                subp.parse(varargin{:});
                sub_pout = subp.Results;
                
                no_days_of_week = sub_pout.no_dow;
                no_time_period = sub_pout.no_time_period;
                no_product = sub_pout.no_product;
                
                % first the time periods
                l = gobjects(0,1);
                legend_cell = {};
                if ~no_time_period
                    fns = fieldnames(time_period_colors);
                    for i_fn = 1:numel(fns)
                        this_tp = time_period_colors.(fns{i_fn});
                        if this_tp.used
                            l(end+1) = line(nan, nan, 'linewidth', connector_linewidth, 'color', this_tp.color, 'parent', parent);
                            legend_cell{end+1} = this_tp.name;
                        end
                    end
                end
                
                if ~no_days_of_week
                    fns = fieldnames(dow_markers);
                    for i_fn = 1:numel(fns)
                        this_dow = dow_markers.(fns{i_fn});
                        if this_dow.used
                            l(end+1) = line(nan,nan,'marker', this_dow.marker, 'linestyle', 'none', 'color', 'k', 'linewidth', marker_linewidth, 'markersize', marker_size, 'parent', parent);
                            legend_cell{end+1} = this_dow.name;
                        end
                    end
                end
                
                if ~no_product
                    fns = fieldnames(product_series);
                    which_products_used = nan(1, numel(fns));
                    for i_fn = 1:numel(fns)
                        which_products_used(i_fn) = product_series.(fns{i_fn}).used;
                    end
                    % We only want to indicate the different products in the
                    % legend if more than one was used, otherwise it is
                    % extraneous information
                    if sum(which_products_used) > 1
                        for i_fn = 1:numel(fns)
                            this_product = product_series.(fns{i_fn});
                            if this_product.used
                                l(end+1) = line(nan,nan,this_product.style);
                                legend_cell{end+1} = this_product.name;
                            end
                        end
                    end
                end
                
                if ~isempty(legend_cell)
                    legend(parent, l', legend_cell, 'location', 'best'); 
                end
            end
            
            function fmt = make_connector_fmt(group_fmt, is_change_sig, are_fits_good)
                these_time_period_colors = {group_fmt.color};
                
                if is_change_sig
                    linestyle = '-';
                else
                    linestyle = '--';
                end 
                
                if all(cellfun(@(x) isequal(x, these_time_period_colors{1}), these_time_period_colors))
                    connector_color = these_time_period_colors{1};
                else
                    connector_color = 'k';
                    % Quick kludge to avoid connecting all days-of-week
                    % points that don't have the same pairwise relationship
                    % that weekend-weekdays do. Since color == time period,
                    % if the colors differ, we don't connect.
                    linestyle = 'none';
                end
                
                % If one or both fits are bad, we should hide this line if
                % the plotting is set to not display bad fits, otherwise
                % leave it, since the point will be plotted in gray.
                if any(~are_fits_good)
                    if strcmpi(bad_fit_display, 'no')
                        linestyle = 'none';
                    else
                        linestyle = ':';
                    end
                end
                
                fmt = struct('color', connector_color, 'linestyle', linestyle, 'linewidth', connector_linewidth);
            end
            
            function group_fmt = make_group_fmt(group_fmt, are_fits_good)
                if strcmpi(bad_fit_display, 'no')
                    fmt_field = 'marker';
                    fmt_val = 'none';
                elseif any(strcmpi(bad_fit_display, {'grey','gray'}))
                    fmt_field = 'color';
                    fmt_val = [0.8 0.8 0.8];
                end
                
                if ~are_fits_good(1)
                    group_fmt(1).(fmt_field) = fmt_val;
                end
                if ~are_fits_good(2)
                    group_fmt(2).(fmt_field) = fmt_val;
                end
            end
            %     %
            % End %
            %     %
        end
        
        function [figs, info] = plot_vcd_vs_concentration(varargin)
            % This plots the surface or boundary layer concentration
            % inferred from average WRF-Chem profiles vs. BEHR VCDs. We
            % will assume that the WRF-Chem profiles are representative of
            % the all days-of-week average VCDs (since the NEI emissions do
            % not have a weekend effect, I'm not sure at the moment if they
            % represent an average of all days or are a representative
            % weekday inventory). So to get at the TWRF and US VCD ->
            % concentrations, we'll interpolate the NO2 profiles by their
            % corresponding VCDs.
            
            p = advInputParser;
            p.addParameter('loc_indices', 1:71);
            p.addParameter('interp_method', 'linear');
            p.addParameter('conc_type', 'boundary layer');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_ind = pout.loc_indices;
            interp_method = pout.interp_method;
            conc_type = pout.conc_type;
            
            locs = misc_emissions_analysis.read_locs_file;
            % Since we're reading the locations directly from the
            % spreadsheet, rather than a .mat file, we can cut it down
            % directly.
            locs = locs(loc_ind);
            
            time_periods = {'beg_2yr','beginning','end_2yr','end'};
            days_of_week = {'UMTWRFS', 'TWRF', 'US'};
            
            n_time = numel(time_periods);
            n_dow = numel(days_of_week);
            n_locs = numel(locs);
            n_levels = [];
            
            vcds = cell(n_time, n_dow);
            profs = cell(n_time, n_dow);
            pres = cell(n_time, n_dow);
            
            for i_time = 1:n_time
                [profs{i_time, 1}, pres{i_time, 1}] = misc_emissions_analysis.avg_wrf_prof_around_loc(locs, time_periods{i_time}, 'nox_or_no2', 'nox', 'radius', 0.35);
                for i_dow = 1:n_dow
                    vcds{i_time, i_dow} = misc_emissions_analysis.avg_vcds_around_loc(locs, time_periods{i_time}, days_of_week{i_dow}, 'radius', 0.35);
                    if i_dow > 1
                        profs{i_time, i_dow} = nan(size(profs{i_time, 1}));
                        pres{i_time, i_dow} = nan(size(pres{i_time, 1}));
                    end
                end
                
                if isempty(n_levels)
                    n_levels = size(profs{i_time, 1},1);
                else
                    if size(profs{i_time, 1}, 1) ~= n_levels
                        error('Different numbers of levels')
                    end
                end
            end
            
            % Now we handle the interpolation and at the same time
            % calculate the concentration, either averaged over the
            % boundary layer or right at the surface
            base_dow_ind = 1;
            concentrations = cell(size(vcds));
            concentrations(:) = {nan(1, n_locs)};
            figs = gobjects(n_locs,1);
            
            dow_markers = misc_emissions_analysis.dow_markers;
            time_period_colors = misc_emissions_analysis.time_period_colors;
            
            for i_loc = 1:n_locs
                figs(i_loc) = figure;
                loc_vcds = cellfun(@(v) v(i_loc), vcds);
                loc_base_sfs = nan(n_levels, n_time);
                %loc_base_pres = nan(n_levels, n_time);
                for i_time = 1:n_time
                    loc_base_sfs(:, i_time) = profs{i_time, base_dow_ind}(:, i_loc);
                    %loc_base_pres(:, i_time) = pres{i_time, base_dow_ind}(:, i_loc);
                end
                
                % For the days-of-week that we haven't declared equivalent
                % to WRF, we need to interpolate to get their shape factors
                for i_time = 1:n_time
                    for i_dow = 1:n_dow
                        if i_dow == base_dow_ind
                            this_prof = profs{i_time, i_dow}(:, i_loc);
                            %this_pres = pres{i_time, i_dow}(:, i_loc);
                        else
                            this_prof = interp_profiles(loc_vcds(:, base_dow_ind), loc_base_sfs, loc_vcds(i_time, i_dow), interp_method);
                            %this_pres = interp_profiles(loc_vcds(:, base_dow_ind), loc_base_pres, loc_vcds(i_time, i_dow), 'nearest');
                        end
                        
                        nox_conc_prof = this_prof .* loc_vcds(i_time, i_dow);
                        if strcmpi(conc_type, 'boundary layer')
                            concentrations{i_time, i_dow}(i_loc) = avg_bl_conc(nox_conc_prof);
                        elseif strcmpi(conc_type, 'surface')
                            concentrations{i_time, i_dow}(i_loc) = nox_conc_prof(1);
                        end
                        
                        marker = dow_markers.(days_of_week{i_dow}).marker;
                        color = time_period_colors.(time_periods{i_time}).color;
                        line(loc_vcds(i_time, i_dow), concentrations{i_time, i_dow}(i_loc), 'color', color, 'marker', marker, 'linestyle', 'none', 'linewidth', 2, 'markersize', 16);
                    end
                end
                
                info.concentrations = concentrations;
                info.vcds = vcds;
                info.locs = locs;
                
            end
            
            function interp_profs = interp_profiles(base_vcds, base_profs, target_vcds, method)
                interp_profs = nan(size(base_profs,1), 1);
                for i_lev = 1:numel(interp_profs)
                    interp_profs(i_lev) = interp1(base_vcds, base_profs(i_lev, :), target_vcds, method, 'extrap');
                end
            end
            
            function conc = avg_bl_conc(nox_prof)
                % Assume the first ~10 model levels are in the BL, based on
                % Chicago
                conc = nanmean(nox_prof(1:10));
            end
        end
        
        function plot_wrf_emissions(varargin)
            E = JLLErrors; 
            p = inputParser;
            p.addParameter('locs',nan); % default value of nan instead of empty because often I use empty to mean do not filter locations, but here nan means I need to ask about that
            p.addParameter('years',[]);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            locs = pout.locs;
            emis_years = pout.years;
            
            if isnan(locs)
                if ask_yn('Plot around specific locations?')
                    locs = misc_emissions_analysis.get_loc_inds_interactive();
                else
                    locs = [];
                end
            end
            
            if isnumeric(locs)
                loc_inds = locs;
                locs = misc_emissions_analysis.read_locs_file();
                locs = locs(loc_inds);
            elseif ~isstruct(locs)
                E.badinput('"loc" must be a scalar number or structure');
            end
            
            if isempty(emis_years)
                emis_years = ask_number('Enter the years to map (separated by space if more than 1)', 'testfxn', @misc_emissions_analysis.is_year_valid, 'testmsg', 'Only years between 2005 and 2014 valid');
            elseif ~isnumeric(emis_years) || ~misc_emissions_analysis.is_year_valid(emis_years)
                [~, valid_years_str] = misc_emissions_analysis.is_year_valid([]);
                E.badinput('"emis_years" must be numeric and in years %s', valid_years_str);
            end
            
            [avg_nei_emissions, nei_lon, nei_lat] = misc_emissions_analysis.load_nei_by_year(emis_years);
            
            if isempty(locs)
                figure;
                pcolor(nei_lon, nei_lat, avg_nei_emissions);
                shading flat
                cb = colorbar;
                cb.Label.String = 'NEI NO Emissions (Mg NO/hr)';
                state_outlines('w');
            else
                % Need to convert the radius to degrees. Assume ~110 km per
                % degree
                km2deg = 1/110;
                for i_loc = 1:numel(locs)
                    [xx,yy] = misc_emissions_analysis.find_indices_box_by_loc_radius(locs(i_loc), locs(i_loc).Radius * km2deg * 3, nei_lon, nei_lat);
                    figure;
                    pcolor(nei_lon(xx,yy), nei_lat(xx,yy), avg_nei_emissions(xx,yy));
                    line(locs(i_loc).Longitude, locs(i_loc).Latitude, 'color', 'w', 'marker', 'p');
                    title(locs(i_loc).Location);
                    cb = colorbar;
                    cb.Label.String = 'NEI NO Emissions (Mg NO/hr)';
                    state_outlines('w');
                end
            end
            
        end
        
        function plot_wind_distribution(varargin)
            p = inputParser;
            p.addParameter('winds_file','');
            p.addParameter('all_winds', nan);
            p.addParameter('loc_inds',nan);
            p.addParameter('as_fraction',nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            use_all_winds = pout.all_winds;
            winds_file = pout.winds_file;
            loc_inds = pout.loc_inds;
            as_fraction = pout.as_fraction;
            
            if isempty(winds_file)
                winds_file = misc_emissions_analysis.get_line_dens_file_interactive();
            end
            
            if isnan(use_all_winds)
                use_all_winds = ask_yn('Use all winds (not just those actually used by calc_line_density)?');
            end
            
            if isnan(loc_inds)
                loc_inds = misc_emissions_analysis.get_loc_inds_interactive();
            end
            
            if isnan(as_fraction)
                as_fraction = ask_yn('Plot bins as fraction (rather than counts)?');
            end
            
            speed_cutoff = 3;
            locs_slow = misc_emissions_analysis.bin_wind_distribution(winds_file, 'all_winds', use_all_winds, 'loc_inds', loc_inds, 'wind_op','lt','wind_speed',speed_cutoff);
            locs_fast = misc_emissions_analysis.bin_wind_distribution(winds_file, 'all_winds', use_all_winds, 'loc_inds', loc_inds, 'wind_op','gt','wind_speed',speed_cutoff);
            
            for i_loc = 1:numel(locs_slow)
                y_val_fast = locs_fast(i_loc).WindDirBinCounts;
                y_val_slow = locs_slow(i_loc).WindDirBinCounts;
                x_val = locs_fast(i_loc).WindDirBinCenters; % this should be the same in both
                y_label_string = '# orbits';
                if as_fraction
                    y_val_fast = y_val_fast ./ sum(locs_fast(i_loc).WindDirBinCounts);
                    y_val_slow = y_val_slow ./ sum(locs_slow(i_loc).WindDirBinCounts);
                    y_label_string = 'Fraction of orbits';
                end
                
                [~,~,~,fit_data] = calc_fit_line(y_val_slow, y_val_fast, 'regression', 'RMA');
                
                figure;
                plot(x_val, y_val_slow, 'bo', x_val, y_val_fast, 'rx');
                legend(sprintf('Wind speed < %d', speed_cutoff), sprintf('Wind speed > %d', speed_cutoff));
                xlabel('Wind direction (degrees CCW of east)');
                ylabel(y_label_string);
                title(locs_fast(i_loc).Location);
                x_limits = get(gca,'xlim');
                y_limits = get(gca,'ylim');
                text(interp1([0 1], x_limits, 0.8), interp1([0 1], y_limits, 0.8), sprintf('R^2 = %.2f', fit_data.R2), 'fontsize', 16);
            end
            
        end
        
        function [wind_bin_weights, wind_bin_edges] = calculate_wind_bin_weights(locs, wind_speed, varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('all_winds', false);
            p.addParameter('bin_width',45);
            p.parse(varargin{:});
            pout = p.Results;
            
            use_all_winds = pout.all_winds;
            bin_width = pout.bin_width;
            
            wind_bin_weights = cell(size(locs));
            wind_bin_edges = cell(size(locs));
            for i_loc = 1:numel(locs)
                % include 'loc_inds' = [] to indicate that we do not want
                % to cut down the locations at all.
                locs_tmp = misc_emissions_analysis.bin_wind_distribution(locs(i_loc), 'all_winds', use_all_winds, 'loc_inds', [], 'wind_op', 'lt', 'wind_speed', wind_speed, 'bin_width', bin_width);
                slow_wind_counts = locs_tmp.WindDirBinCounts;
                locs_tmp = misc_emissions_analysis.bin_wind_distribution(locs(i_loc), 'all_winds', use_all_winds, 'loc_inds', [], 'wind_op', 'gt', 'wind_speed', wind_speed, 'bin_width', bin_width);
                fast_wind_counts = locs_tmp.WindDirBinCounts;
                wind_bin_edges{i_loc} = locs_tmp.WindDirBinEdges;
                
                wind_bin_weights{i_loc} = fast_wind_counts ./ slow_wind_counts;
            end
            
            
        end
        
        function locs = bin_wind_distribution(winds_file, varargin)
            % BIN_WIND_DISTRIBUTION Bin the occurrences of wind directions
            %   BIN_WIND_DISTRIBUTION( WINDS_FILE ) Given a file name as a
            %   char array, will load the file WINDS_FILE and try to bin
            %   the data within it. The file must contain a structure
            %   called "locs" with fields "WindDir" and "WindSpeed".
            %   Alternately, give the locs structure contained in such a
            %   file as WINDS_FILE directly. Returns the locs structure
            %   with the fields "WindDirBinCounts", "WindDirBinEdges", and
            %   "WindDirBinCenters" fields added.
            %
            %   Parameters:
            %       'all_winds' - boolean, whether to include only winds
            %       for orbits used by calc_line_density() (false, default)
            %       or all winds (true). If false, then the locs structure
            %       in the winds file must also contain the field
            %       "WindUsedBool"
            %
            %       'loc_inds' - a logical or numeric index array to cut
            %       down the locs structure to only relevant locations.
            %       Asked interactively if omitted (or given as NaN).
            %
            %       'wind_op' - either the string 'lt' or 'gt' to indicate
            %       whether only winds less than or greater than the
            %       criterion speed should be used. Asked interactively if
            %       omitted.
            %
            %       'wind_speed' - the wind speed criterion (in
            %       meters/second) that goes along with 'wind_op'. Asked
            %       interactively if omitted.
            %
            %       'bin_width' - the width of the bins in degrees. Default
            %       is 45.
            E = JLLErrors;
            p = inputParser;
            p.addParameter('all_winds', false);
            p.addParameter('loc_inds',nan);
            p.addParameter('wind_op','');
            p.addParameter('wind_speed',[]);
            p.addParameter('bin_width',45);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            use_all_winds = pout.all_winds;
            loc_inds = pout.loc_inds;
            wind_op = pout.wind_op;
            wind_speed = pout.wind_speed;
            bin_width = pout.bin_width;
            
            if ischar(winds_file)
                W = load(winds_file);
                locs = W.locs;
            elseif isstruct(winds_file)
                locs = winds_file;
            end
            if any(~isfield(locs, {'WindSpeed','WindDir'}))
                E.badinput('The file "%s" does not appear to have the expected wind data', winds_file);
            end
            
            if ~use_all_winds && ~isfield(locs, 'WindUsedBool')
                if ischar(winds_file)
                    msg = sprintf('The file "%s" does not contain the field "WindUsedBool" in the locs structure. Pass "''all_wind'', true" if you need to bin a file that does not have this field', winds_file);
                else
                    msg = 'The given structure does not contain the field "WindUsedBool". Pass "''all_wind'', true" if you need to bin a file that does not have this field';
                end
                E.callError('missing_wind_used', msg)
            end
            
            if isnan(loc_inds)
                loc_inds = misc_emissions_analysis.get_loc_inds_interactive();
            end
            if ~isempty(loc_inds)
                locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            end
            
            if ~isnumeric(bin_width) || ~isscalar(bin_width) || bin_width <= 0
                E.badinput('''bin_width'' must be a positive scalar number');
            end
            
            [wind_op, wind_speed] = misc_emissions_analysis.choose_wind_criteria(wind_op, wind_speed);
            
            % TODO: find a way to extract winds actually used in the
            % analysis, i.e. those that match up with BEHR swaths that have
            % actual data
            for i_loc = 1:numel(locs)
                winds_logical = misc_emissions_analysis.set_wind_conditions(locs(i_loc),wind_speed,wind_op,'none');
                if ~use_all_winds
                    winds_logical = winds_logical & locs(i_loc).WindUsedBool;
                end
                all_wind_directions = locs(i_loc).WindDir(winds_logical);
                [bin_counts, bin_edges] = histcounts(all_wind_directions, -180:bin_width:180);
                locs(i_loc).WindDirBinCounts = bin_counts;
                locs(i_loc).WindDirBinEdges = bin_edges;
                locs(i_loc).WindDirBinCenters = 0.5*(bin_edges(1:end-1)+bin_edges(2:end));
            end
            
        end
        
        function plot_loc_wind_rose(varargin)
            E=JLLErrors;
            p = inputParser;
            p.addParameter('time_period', '');
            p.addParameter('loc_inds', nan);
            p.addParameter('wind_op','');
            p.addParameter('wind_speed',[]);
            p.addParameter('convert_wind_def', true);
            p.addParameter('keep_rejected_directions', nan);
            p.addParameter('use_wrf', nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            time_period = pout.time_period;
            loc_inds = pout.loc_inds;
            wind_op = pout.wind_op;
            wind_speed = pout.wind_speed;
            convert_wind_def_bool = pout.convert_wind_def;
            keep_rejected_directions_bool = pout.keep_rejected_directions;
            use_wrf = pout.use_wrf;
            
            [start_dates, end_dates] = misc_emissions_analysis.select_start_end_dates(time_period);
            winds_file = misc_emissions_analysis.winds_file_name(start_dates{1}, end_dates{end});
            W = load(winds_file);
            W.locs = misc_emissions_analysis.append_new_spreadsheet_fields(W.locs);
            
            if isnan(loc_inds)
                loc_inds = misc_emissions_analysis.get_loc_inds_interactive();
            end
            if ~isempty(loc_inds)
                W.locs = misc_emissions_analysis.cutdown_locs_by_index(W.locs, loc_inds);
            end
            
            [wind_op, wind_speed] = misc_emissions_analysis.choose_wind_criteria(wind_op, wind_speed);
            
            if isnan(keep_rejected_directions_bool)
                keep_rejected_directions_bool = ask_yn('Retain rejected wind directions for the wind rose?');
            end
            
            if ~keep_rejected_directions_bool
                [~, ~, wind_rej_field] = misc_emissions_analysis.ask_to_use_wrf(use_wrf);
            else
                wind_rej_field = 'none';
            end
            
            
            for i_loc = 1:numel(W.locs)
                winds_logical = misc_emissions_analysis.set_wind_conditions(W.locs(i_loc), wind_speed, wind_op, wind_rej_field);
                wind_dir = W.locs.WindDir(winds_logical);
                wind_vel = W.locs.WindSpeed(winds_logical);
                if convert_wind_def_bool
                    % In my scheme, a wind direction of 0 degrees means the
                    % wind is blowing to the east, 90 means to the north,
                    % 180 to the west, and -90 to the south. Since a wind
                    % rose traditionally represents the direction that the
                    % wind is coming from, we need to flip this. That means
                    % a wind coming from the north is blowing to the south,
                    % so in my definition a northerly wind has an angle of
                    % -90. An easterly has an angle of +/- 180 in my
                    % definition.
                    north_angle = -90;
                    east_angle = -180;
                else
                    north_angle = 90;
                    east_angle = 0;
                    fprintf('!!! NOTE: wind rose plotting directions winds are blowing TOWARDS !!!\n');
                end
                
                % Wind rose automatically creates a new figure
                WindRose(wind_dir, wind_vel, 'anglenorth', north_angle, 'angleeast', east_angle);
            end
        end
        
        function save_tables()
            time_periods = {'beg_2yr', 'beginning', 'end_2yr', 'end'};
            req_values = {'FitAccepted', 'Correlation', 'FracPtsBiasedWindow10', 'FracPtsBiasedWindow20', 'LifetimesDownwind', 'R2', 'Tau', 'TauPercentUncertainty', 'E', 'EPercentUncertainty'};
            
            common_opts = {'time_periods', time_periods, 'values', req_values};
            BEHRTables = misc_emissions_analysis.tabulate_parameter_correlation(common_opts{:}, 'days_of_week', {'UMTWRFS','TWRF','US'}, 'data', 'BEHR');
            WRFTables = misc_emissions_analysis.tabulate_parameter_correlation(common_opts{:}, 'days_of_week', {'UMTWRFS'}, 'data', 'BEHR');
            savename = fullfile(misc_emissions_analysis.table_dir, sprintf('Tables%s.mat', datestr(today, 'yyyymmdd')));
            save(savename, 'BEHRTables', 'WRFTables');
        end
        
        function output_tables = tabulate_parameter_correlation(varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('time_periods', {});
            p.addParameter('days_of_week', {});
            p.addParameter('data', '')
            p.addParameter('values', {});
            
            p.parse(varargin{:});
            pout = p.Results;
            
            quantity_fxns = struct('Correlation', @calc_correlation,...
                'FitAccepted', @(this_loc) misc_emissions_analysis.is_fit_good_by_loc(this_loc),...
                'R2', @(this_loc) this_loc.fit_info.param_stats.r2,...
                'FracPtsBiasedWindow10', @(this_loc) mean(misc_emissions_analysis.test_fit_for_bias(this_loc.no2_sectors.linedens, this_loc.fit_info.emgfit, 'window', 10)),...
                'FracPtsBiasedWindow20', @(this_loc) mean(misc_emissions_analysis.test_fit_for_bias(this_loc.no2_sectors.linedens, this_loc.fit_info.emgfit, 'window', 20)),...
                'LifetimesDownwind', @(this_loc) misc_emissions_analysis.n_lifetimes_downwind(this_loc.no2_sectors.x, this_loc.fit_info.ffit.x_0, this_loc.fit_info.ffit.mu_x),...
                'E', @(this_loc) this_loc.emis_tau.emis,...
                'EPercentUncertainty', @(this_loc) this_loc.emis_tau.emis_uncert / this_loc.emis_tau.emis * 100,...
                'Tau', @(this_loc) this_loc.emis_tau.tau,...
                'TauPercentUncertainty', @(this_loc) this_loc.emis_tau.tau_uncert / this_loc.emis_tau.tau * 100);
            
            avail_time_periods = {'beg_2yr', 'beginning', 'end_2yr', 'end'};
            time_periods = opt_ask_multiselect('Which time periods to include?', avail_time_periods, pout.time_periods, '"time_periods"');
            days_of_week = opt_ask_multiselect('Which days-of-week to include?', {'UMTWRFS','TWRF','US'}, pout.days_of_week, '"days_of_week"');
            data_source = opt_ask_multichoice('Which data set to use?', {'BEHR', 'WRF'}, pout.data, '"data"', 'list', true);
            requested_values = opt_ask_multiselect('Which values to tabulate?', fieldnames(quantity_fxns), pout.values, '"values"');
            requested_values = cellfun(@capitalize_words, requested_values, 'uniform', false);
            requested_values = regexprep(requested_values, '\s', '');
            
            wrf_bool = strcmpi(data_source, 'WRF');
            ss_locs = misc_emissions_analysis.read_locs_file('Cities','PowerPlants');
            ss_loc_names = {ss_locs.ShortName};
            n_locs = numel(ss_locs);
            
            n_columns = numel(time_periods)*numel(days_of_week);
            values = nan(n_locs, n_columns);
            colnames = cell(1, n_columns);
            
            output_tables = make_empty_struct_from_cell(requested_values, values);
            
            
            i_column = 1;
            for i_time = 1:numel(time_periods)
                fprintf('Working on %s\n', time_periods{i_time});
                if regcmp(time_periods{i_time}, '2yr') && ~wrf_bool
                    loc_inds = 1:70;
                else
                    loc_inds = 1:71;
                end
                
                for i_dow = 1:numel(days_of_week)
                    fprintf('  Loading %s data\n', days_of_week{i_dow});
                    
                    [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_periods{i_time});
                    Data = load(misc_emissions_analysis.fits_file_name(start_date, end_date, wrf_bool, loc_inds, days_of_week{i_dow}, 'lu'));
                    
                    for i_val = 1:numel(requested_values)
                        fprintf('  Computing %s\n', requested_values{i_val});
                        for i_loc = 1:numel(Data.locs)
                            this_loc = Data.locs(i_loc);
                            if isempty(this_loc.fit_info.emgfit)
                                fprintf('    Skipping %s because it looks like the fit didn''t work\n', this_loc.ShortName);
                                continue
                            end
                            loc_ind = strcmp(this_loc.ShortName, ss_loc_names);
                            
                            output_tables.(requested_values{i_val})(loc_ind, i_column) = quantity_fxns.(requested_values{i_val})(this_loc);
                        end
                    end
                    
                    colnames{i_column} = sprintf('%s_%s', time_periods{i_time}, days_of_week{i_dow});
                    i_column = i_column + 1;
                end
            end
            
            for i_val = 1:numel(requested_values)
                output_tables.(requested_values{i_val}) = array2table(output_tables.(requested_values{i_val}), 'RowNames', ss_loc_names, 'VariableNames', colnames);
            end
            
            function corr_val = calc_correlation(this_loc)
                % for now, just want correlation of a and x_0
                var1_ind = 1;
                var2_ind = 2;
                try
                    [~, corr_mat] = cov2corr(emg_cov_mat(this_loc.fit_info.fitresults.Hessian, this_loc.no2_sectors.x, this_loc.no2_sectors.linedens, this_loc.fit_info.emgfit));
                    corr_val =  corr_mat(var1_ind, var2_ind);
                catch err
                    if strcmpi(err.identifier, 'finance:cov2corr:invalidCovMatrixSymmetry')
                        fprintf('    Skipping %s due to problem calculating correlation matrix: "%s"\n', this_loc.ShortName, err.message);
                        corr_val = nan;
                    else
                        rethrow(err)
                    end
                end
            end
        end
        
        
        function bias_found = test_fit_for_bias(line_dens, emgfit, varargin)
            p = advInputParser;
            p.addParameter('window', 10);
            p.addParameter('confidence', 0.95);
            
            p.parse(varargin{:})
            pout = p.Results;
            
            window = pout.window;
            confidence = pout.confidence;
            
            % tinv assumes a one-sided distribution, the two sided one has
            % half the alpha (1-p) value
            alpha = 1 - confidence;
            p = 1 - alpha/2;
            t_value = tinv(p, window);
            
            bias_found = false(size(line_dens));
            
            for i_start = 1:(numel(line_dens) - window + 1)
                i_end = i_start + window - 1;
                mean_diff = abs(nanmean(emgfit(i_start:i_end)) - nanmean(line_dens(i_start:i_end)));
                diff_crit = t_value * nanstd(line_dens(i_start:i_end)) / sqrt(window);
                if mean_diff > diff_crit
                    bias_found(i_start:i_end) = true;
                end
            end
        end
    end
    
    methods(Static = true, Access = private)
        function [bool, years_list] = is_year_valid(years_in)
            bool = all(years_in == 2005 | (years_in >= 2007 & years_in <= 2009) | (years_in >= 2012 & years_in <= 2014));
            years_list = '2005, 2007-09, 2012-14';
        end
        
        function [first_time_period, second_time_period, first_weekdays, second_weekdays, loc_types, series_labels] = get_change_file_selection_input(varargin)
            allowed_change_types = {'decadal','weekend-weekday'};
            allowed_time_periods = {'beginning','end'};
            allowed_loc_types = {'Cities','PowerPlants'};
            
            p = inputParser;
            p.addParameter('change_type', '', @(x) ischar(x) && ismember(x, [{''}, allowed_change_types]));
            p.addParameter('time_period', '', @(x) ischar(x) && ismember(x, [{''}, allowed_time_periods]));
            p.addParameter('loc_types','');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            change_type = pout.change_type;
            time_period = pout.time_period;
            loc_types = pout.loc_types;
            
            if isempty(change_type)
                change_type = ask_multichoice('Which type of change to plot?', allowed_change_types, 'list', true);
            end
            
            if ~strcmpi(change_type, 'decadal')
                if isempty(time_period)
                    time_period = ask_multichoice('For which time period?', allowed_time_periods, 'list', true);
                end
                
                first_time_period = time_period;
                second_time_period = time_period;
                first_weekdays = 'TWRF';
                second_weekdays = 'US';
                series_labels = {'Weekday','Weekend'};
            else
                first_time_period = 'beginning';
                second_time_period = 'end';
                first_weekdays = 'UMTWRFS';
                second_weekdays = 'UMTWRFS';
                series_labels = {'2005-2007','2012-2013'};
            end
            
            if ~ischar(loc_types) && ~iscellstr(loc_types)
                E.badinput('"loc_types" must be a character array or cell array of character arrays')
            elseif ~isempty(loc_types)
                if ~all(ismember(loc_types, allowed_loc_types))
                    E.badinput('"loc_types" must be one or more of the following: %s', strjoin(allowed_loc_types, ', '));
                end
            else
                loc_types = ask_multiselect('Select one or more site types to plot:', allowed_loc_types);
            end
        end
        
        function avg_vcds = avg_vcds_around_loc(locs, time_period, days_of_week, varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('radius', 'by_loc');
            p.addParameter('species', 'no2');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_radius = pout.radius;
            vcd_species = pout.species;
            
            if ischar(avg_radius)
                if strcmpi(avg_radius, 'by_loc')
                    avg_radius = [];
                else
                    E.badinput('The only valid string for "avg_radius" is "by_loc"')
                end
            end
            
            allowed_species = {'no2','hcho'};
            if ~ismember(vcd_species, allowed_species)
                E.badinput('"vcd_species" must be one of: %s', strjoin(allowed_species, ', '))
            elseif strcmp(vcd_species, 'hcho')
                % currently I only have all day-of-week HCHO VCD averages
                % because I don't expect large weekday-weekend difference
                % (but might be worth checking) so for now we'll just force
                % the HCHO VCDs to be from all days of week
                
                if ~strcmpi(days_of_week, 'UMTWRFS')
                    warning('HCHO columns will be from all days of week, even if a different subset was specified');
                end
                days_of_week = 'UMTWRFS';
            end
            
            [start_dates, end_dates] = misc_emissions_analysis.select_start_end_dates(time_period);
            time_period_years = unique(cellfun(@year, veccat(start_dates, end_dates)));
            VCDs = load(misc_emissions_analysis.avg_file_name(time_period_years, days_of_week, vcd_species));
            lon_res = mean(diff(VCDs.daily.lon(1,:)));
            lat_res = mean(diff(VCDs.daily.lat(:,1)));
            if abs(lon_res - lat_res) > 1e-10
                E.notimplemented('Different lon and lat resolutions');
            end
            avg_vcds = nan(size(locs));
            for i_loc = 1:numel(locs)
                % This will use the box width (from center to edge
                % perpendicular to the wind direction) as the radius and
                % find all grid points with centers within that radius.
                xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(locs(i_loc), VCDs.daily.lon, VCDs.daily.lat, avg_radius);
                avg_vcds(i_loc) = nanmean(VCDs.daily.(vcd_species)(xx_radius));
            end
        end
        
        function [shape_factors, pres_levels] = avg_wrf_prof_around_loc(locs, time_period, varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('nox_or_no2', 'nox');
            p.addParameter('radius', 'by_loc');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            nox_or_no2 = pout.nox_or_no2;
            avg_radius = pout.radius;
            
            Profs = misc_emissions_analysis.load_time_averaged_wrf_profs(time_period);
            
            shape_factors = nan(size(Profs.pres,1), numel(locs));
            pres_levels = nan(size(Profs.pres,1), numel(locs));
            
            if ischar(avg_radius)
                if strcmpi(avg_radius, 'by_loc')
                    avg_radius = [];
                else
                    E.badinput('The only valid string for "avg_radius" is "by_loc"')
                end
            end
            
            for i_loc = 1:numel(locs)
                xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(locs(i_loc), Profs.lon, Profs.lat, avg_radius);
                this_no2 = nanmean(Profs.no2(:,xx_radius),2) * 1e-6;
                this_no = nanmean(Profs.no(:,xx_radius),2) * 1e-6;
                this_pres = nanmean(Profs.pres(:,xx_radius),2);
                
                
                % Now we need to calculate the NO2 VCDs, then either the
                % NOx or NO2 shape factors
                no2_vcd = integPr2(this_no2, this_pres, this_pres(1), this_pres(end));
                if strcmpi(nox_or_no2, 'nox')
                    shape_factors(:, i_loc) = (this_no + this_no2) ./ no2_vcd;
                elseif strcmpi(nox_or_no2, 'no2')
                    shape_factors(:, i_loc) = this_no2 ./ no2_vcd;
                else
                    E.notimplemented('No method set for "nox_or_no2" == "%s"', nox_or_no2);
                end
                pres_levels(:, i_loc) = this_pres;
            end
        end
        
        function VCDs = avg_wrf_vcds_around_loc(locs, time_period, specie, varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('radius', 'by_loc');
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_radius = pout.radius;
            
            if ~ischar(specie)
                E.badinput('SPECIE must be a char array')
            end
            
            Profs = misc_emissions_analysis.load_time_averaged_wrf_profs(time_period);
            xx_bad_species = ~ismember(specie, fieldnames(Profs));
            if any(xx_bad_species)
                E.badinput('The following species are not available in the time-averaged profile file: %s', strjoin(specie(xx_bad_species), ', '))
            end
            
            if ischar(avg_radius)
                if strcmpi(avg_radius, 'by_loc')
                    avg_radius = [];
                else
                    E.badinput('The only valid string for "avg_radius" is "by_loc"')
                end
            end
            
            VCDs = nan(size(locs));
            
            for i_loc = 1:numel(locs)
                xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(locs(i_loc), Profs.lon, Profs.lat, avg_radius);
                this_pres = Profs.pres(:,xx_radius);
                this_profiles = Profs.(specie)(:,xx_radius);
                these_vcds = nan(1, size(this_profiles, 2));
                for i_prof = 1:numel(these_vcds)
                    % The species concentrations in WRF are usually in
                    % ppm and the time averaging doesn't change that.
                    % integPr2 requires the straight mixing ratio
                    % (parts-per-part)
                    these_vcds(i_prof) = integPr2(this_profiles(:,i_prof)*1e-6, this_pres(:,i_prof), this_pres(1,i_prof), this_pres(end,i_prof));
                end
                VCDs(i_loc) = nanmean(these_vcds);
            end
        end
        
        
        function Profs = load_time_averaged_wrf_profs(time_period)
            [start_dates, end_dates] = misc_emissions_analysis.select_start_end_dates(time_period);
            prof_file = misc_emissions_analysis.wrf_avg_prof_file(start_dates, end_dates);
            Profs = load(prof_file);
            
            dvec = make_datevec(start_dates, end_dates);
            % My first run of average profiles added already weighted data
            % to a weighted average when combining the different workers,
            % so the weights were double counted. The weights in question
            % were the number of days averaged on each worker (I used 2
            % workers) so we need to divide by that to bring things back
            % into line.
            weight_correction = 1;%/(numel(dvec)/2);
            %warning('2018-07-16: Correcting average profiles for double-counting weights. For average profiles recalculated after 16 Jul 2018, this correction must be turned off');
            
            
            fns = fieldnames(Profs);
            for i_fn = 1:numel(fns)
                this_fn = fns{i_fn};
                if ismatrix(Profs.(this_fn))
                    % 2D arrays like lat and lon should not be rearranged
                    continue
                elseif ndims(Profs.(this_fn) == 3)
                    % Permute 3D arrays so that the first dimension is
                    % vertical, this will make it easier to subset them.
                    Profs.(this_fn) = permute(Profs.(this_fn), [3 1 2]) * weight_correction;
                else
                    E.notimplemented('Did not expect an array with ndims > 3')
                end
            end
        end
        
        function plot_map_series(ax, map_series, varargin)
            ax.NextPlot = 'add';
            for i=1:numel(map_series)
                % Default scatter size is 36. Double the point size for
                % legibility.
                scatter(ax, map_series(i).lon, map_series(i).lat, 72, map_series(i).change, map_series(i).symbol{:});
            end
            
            
            if ismember(varargin, 'include_names')
                % Put the name further out along the same vector from the
                % center as the point being plotted.
                point_lons = veccat(map_series.lon);
                point_lats = veccat(map_series.lat);
                center_lons = veccat(map_series.center_lon);
                center_lats = veccat(map_series.center_lat);
                dlon = point_lons - center_lons;
                dlat = point_lats - center_lats;
                names = veccat(map_series.names);
                % ensure all given as column vectors
                text(center_lons(:) + 1.5*dlon(:),center_lats(:)+1.5*dlat(:),names(:),'parent',ax);
            end
        end
        
        function difference_is_significant = is_change_significant(values, value_sds, value_dofs)
            % Calculate whether the difference is significant at the 95%
            % confidence level.
            
            difference_is_significant = false(size(values,1),1);
            for i_chng = 1:size(values,1)
                % The EMG fits have 5 fitting parameters.
                if any(isnan(values(i_chng,:))) || any(isnan(value_sds(i_chng,:))) || any(imag(value_sds(i_chng,:)) ~= 0) || any(isnan(value_dofs(i_chng,:)))
                    difference_is_significant(i_chng) = false;
                else
                [~, ~, difference_is_significant(i_chng)] = two_sample_t_test(values(i_chng, 1), value_sds(i_chng, 1).^2 .* value_dofs(i_chng, 1), value_dofs(i_chng, 1) + 5,...
                    values(i_chng, 2), value_sds(i_chng, 2).^2 .* value_dofs(i_chng, 2), value_dofs(i_chng, 2) + 5,...
                    sum(value_dofs(i_chng,:)), 'confidence', 0.95);
                end
            end
        end
        
        function plot_series = create_map_series(values, value_sds, value_dofs, value_r2, coords, corner, loc_names)
            % Helper function to set up the series with the right
            % coordinates and symbols based on their confidence and r2
            % values. "corner" should be a string indicating
            % which corner of the triangle in the plot it should be.
            % Possible values are "center", "top", "bottom-left", and
            % "bottom-right"
            
            E = JLLErrors;
            
            if ~exist('corner', 'var')
                corner = 'center';
            end
            
            % Go ahead and compute the offsets for the corner now
            offset_distance = 0.5; % the triangle's corners will be this distance from the center
            if isnumeric(corner)
                offset = offset_distance * [cosd(corner), sind(corner)];
            else
                switch lower(corner)
                    case 'center'
                        offset = [0 0];
                    case 'top'
                        offset = offset_distance * [cosd(90), sind(90)];
                    case 'bottom-left'
                        offset = offset_distance * [cosd(-150), sind(-150)];
                    case 'bottom-right'
                        offset = offset_distance * [cosd(-30), sind(-30)];
                    case 'bottom'
                        offset = offset_distance * [cosd(-90), sind(-90)];
                    otherwise
                        E.badinput('CORNER must be one of the strings "center", "top", "bottom, "bottom-left", or "bottom-right"');
                end
            end
            
            r2_criterion = 0.9;
            
            xx_bad_r2 = any(value_r2 < r2_criterion, 2);
            
            difference_is_significant = misc_emissions_analysis.is_change_significant(values, value_sds, value_dofs);
            
            xx_not_sig = ~difference_is_significant & ~xx_bad_r2;
            xx_sig = ~xx_not_sig & ~xx_bad_r2;
            
            xx_cell = {xx_sig, xx_not_sig, xx_bad_r2};
            
            default_mat = repmat({[]}, size(xx_cell));
            
            % Use cell arrays of cell arrays so that we can call
            % SCATTER(lon, lat, [], change, symbol{:}) and get filled
            % circles but not filled asterisks or x's. Testing shows that
            % calling 'filled' with * markers makes them just disappear.
            symbols = {{'o','filled'},{'*'},{'x'}};
            
            plot_series = struct('lon', default_mat, 'lat', default_mat, 'change', default_mat, 'percent_change', default_mat, 'symbol', symbols, 'names', {{}});
            for i_series = 1:numel(plot_series)
                plot_series(i_series).center_lon = coords.lon(xx_cell{i_series});
                plot_series(i_series).center_lat = coords.lat(xx_cell{i_series});
                plot_series(i_series).lon = coords.lon(xx_cell{i_series}) + offset(1);
                plot_series(i_series).lat = coords.lat(xx_cell{i_series}) + offset(2);
                
                % Take the difference and percent difference of the second
                % column of values minus the first.
                plot_series(i_series).change = diff(values(xx_cell{i_series},:),1,2);
                % reldiff(A,B) is (A-B)/B
                plot_series(i_series).percent_change = reldiff(values(xx_cell{i_series},2), values(xx_cell{i_series},1)) * 100;
                
                plot_series(i_series).names = loc_names(xx_cell{i_series});
            end
            
        end
        
        function dow = choose_days_of_week(varargin)
            E = JLLErrors;
            
            p = advInputParser;
            p.addOptional('days_of_week','',@ischar);
            p.addFlag('individual');
            
            p.parse(varargin{:});
            pout = p.AdvResults;
            
            days_of_week = pout.days_of_week;
            indiv_bool = pout.individual;
            
            if isempty(days_of_week)
                if indiv_bool
                    allowed_dow = {'S','M','T','W','R','F','S'};
                    dow = ask_multiselect('Choose days of week to include', allowed_dow);
                    dow = strjoin(dow, '');
                else
                    allowed_dow = {'weekdays', 'weekends', 'both'};
                    dow_ans = ask_multichoice('Choose days of week to include', allowed_dow, 'list', true);
                    switch lower(dow_ans)
                        case 'weekdays'
                            dow = 'TWRF';
                        case 'weekends'
                            dow = 'US';
                        case 'both'
                            dow = 'UMTWRFS';
                    end
                end
            elseif ~ischar(days_of_week) || any(~ismember(days_of_week, 'UMTWRFS'))
                E.badinput('DAYS_OF_WEEK must be a character array consisting only of the characters U, M, T, W, R, F, or S');
            else
                dow = days_of_week;
            end
        end
        
        function earth_area = calculate_total_omi_grid_area(lon, lat, resolution)
            E = JLLErrors;
            if numel(lon) ~= numel(lat)
                E.badinput('LON and LAT must have the same number of elements');
            end
            if ~isscalar(resolution) || ~isnumeric(resolution)
                E.badinput('RESOLUTION must be a scalar number');
            end
            
            earth_area = nan(size(lon));
            % The units of VCD are molec/cm^2, so we want area in cm^2.
            % WGS84 is the reference used in GPS systems.
            reference_ellip = referenceEllipsoid('wgs84','cm');
            for i_pt = 1:numel(lon)
                earth_area(i_pt) = areaquad(lat(i_pt)-resolution/2, lon(i_pt)-resolution/2, lat(i_pt)+resolution/2, lon(i_pt)+resolution/2, reference_ellip);
            end
        end
        
        function [xx,yy] = find_indices_box_by_loc_radius(location, radius, lon_grid, lat_grid )
            % Convenience wrapper around find_indicies_in_box_around_point
            % that converts a radius to a number of boxes.
            
            % Need to figure out the average grid spacing around the
            % location. Check each of the four cardinal directions and
            % average them.
            [pt_ind(1), pt_ind(2)] = misc_emissions_analysis.find_lat_lon_index(location.Longitude, location.Latitude, lon_grid, lat_grid);
            pt_moves = [-1 0; 1 0; 0 -1; 0 1];
            grid_del = nan(size(pt_moves,1),1);
            for i_pt = 1:size(pt_moves,1)
                i1 = pt_ind(1);
                j1 = pt_ind(2);
                i2 = pt_ind(1) + pt_moves(1);
                j2 = pt_ind(2) + pt_moves(2);
                grid_del(i_pt) = sqrt((lon_grid(i1,j1) - lon_grid(i2,j2)).^2 + (lat_grid(i1,j1) - lat_grid(i2,j2)).^2);
            end
            
            grid_del = nanmean(grid_del);
            n_cells = ceil(radius / grid_del);
            [xx, yy] = misc_emissions_analysis.find_indicies_in_box_around_point(location, lon_grid, lat_grid, n_cells);
        end
        
        function [xx,yy] = find_lat_lon_index(lon_pt, lat_pt, lon_grid, lat_grid)
            % Finds the index of the grid point closest to lon_pt and
            % lat_pt
            r = sqrt((lon_grid - lon_pt).^2 + (lat_grid - lat_pt).^2);
            [~, i_min] = min(r(:));
            [xx, yy] = ind2sub(size(lon_grid), i_min); 
        end
        
        function locs = cutdown_locs_by_index(locs, loc_inds, varargin)
            % Cuts down the structure LOCS to the locations specified by
            % LOCS_INDS by matching the site names from LOCS against the
            % location names in the structure read in from the spreadsheet.
            % If LOCS_INDS is an empty array, do not cut the locations down
            % at all.
            
            E = JLLErrors;
            p = advInputParser;
            p.addFlag('keep_order');
            p.parse(varargin{:});
            pout = p.Results;
            
            keep_order = pout.keep_order;
            
            if isempty(loc_inds)
                return 
            end
            
            locs_ss = misc_emissions_analysis.read_locs_file();
            ss_names = {locs_ss(loc_inds).Location};
            in_names = {locs.Location};
            if ~keep_order
                xx = ismember(in_names, ss_names);
                if sum(xx) ~= numel(loc_inds)
                    E.callError('loc_lookup_error', 'A different number of locations was found in LOCS that specified by LOC_INDS');
                end
            else
                xx = nan(size(ss_names));
                for i=1:numel(xx)
                    xx_i = find(strcmp(ss_names{i}, in_names));
                    if isempty(xx_i)
                        E.callError('loc_lookup_error', 'Could not find "%s" in the given locs structure', ss_names{i});
                    end
                    xx(i) = xx_i;
                end
            end
            locs = locs(xx);
            
        end
        
        function new_locs = match_locs_structs(new_locs, base_locs)
            % Cut down NEW_LOCS to have the same locations as BASE_LOCS
            E = JLLErrors;
            new_loc_names = {new_locs.Location};
            base_loc_names = {base_locs.Location};
            xx = ismember(new_loc_names, base_loc_names);
            new_locs(~xx) = [];
            check_loc_names = {new_locs.Location};
            if numel(check_loc_names) ~= numel(base_loc_names)
                E.callError('locs_not_available', 'One or more locations in BASE_LOCS were not present in NEW_LOCS');
            elseif any(~strcmp(check_loc_names, base_loc_names))
                E.callError('locs_not_matched', 'Trouble matching the cutdown NEW_LOCS to BASE_LOCS, the locations may be out of order');
            end
        end
        
        function loc_inds = loc_types_to_inds(varargin)
            locs = misc_emissions_analysis.read_locs_file();
            loc_inds = false(size(locs));
            for i_loc = 1:numel(locs)
                loc_inds(i_loc) = ismember(locs(i_loc).SiteType, varargin);
            end
            loc_inds = find(loc_inds);
        end
        
        function locs = append_new_spreadsheet_fields(locs)
            % If I have added new values to the spreadsheet since the last
            % time a particular step was run, this will ensure those values
            % are includes in the locations structure given as input.
            spreadsheet_locs = misc_emissions_analysis.read_locs_file();
            locs = copy_structure_fields(spreadsheet_locs, locs, 'missing');
        end
        
        
        function [is_sectors_vec, is_filtered_vec, is_weighted_vec, winds_strings, is_wrf_vec, days_of_week] = extract_info_from_file_names(files)
            E = JLLErrors;
            if ~ischar(files) && ~iscellstr(files)
                E.badinput('FILES must be a char array or cell array of such');
            end
            is_sectors_vec = regcmp(files, 'sectors');
            is_filtered_vec = ~regcmp(files, 'unfiltered'); % have to use not here b/c regcmp(files, 'filtered') would match unfiltered as well
            is_weighted_vec = ~regcmp(files, 'unweighted');
            winds_strings = regexp(files, 'winds-(lt|gt)\d', 'match', 'once');
            is_wrf_vec = regcmp(files, '^WRF');
            days_of_week = regexp(files, '(?<=_)[UMTWRFS]+(?=\.mat)', 'match', 'once');
        end
        
        function [use_wrf, loc_inds, wind_rej_field, file_loc_inds] = ask_to_use_wrf(use_wrf, default_loc_inds, default_file_inds)
            if nargin < 1 || isnan(use_wrf)
                use_wrf = ask_yn('Use WRF line densities? (BEHR if no)');
            end
            
            if ~exist('default_loc_inds', 'var') || isempty(default_loc_inds)
                default_loc_inds = 1:71;
            end
            if ~exist('default_file_inds', 'var') || isempty(default_file_inds)
                default_file_inds = default_loc_inds;
            end
            
            if use_wrf
                loc_inds = default_loc_inds(ismember(default_loc_inds, 1:71));
                wind_rej_field = misc_emissions_analysis.wind_reject_field_wrf;
                file_loc_inds = 1:71;
            else
                loc_inds = default_loc_inds;
                file_loc_inds = default_file_inds;
                wind_rej_field = misc_emissions_analysis.wind_reject_field_std;
            end
            
        end
        
        function [wind_op, wind_speed] = choose_wind_criteria(wind_op, wind_speed)
            E = JLLErrors;
            
            allowed_wind_ops = {'lt','gt'};
            wind_op_names = {'less than', 'greater than'};
            if isempty(wind_op)
                i_op = ask_multichoice('Use winds that are what vs the criteria?', wind_op_names, 'list', true, 'index', true);
                wind_op = allowed_wind_ops{i_op};
            elseif ~ismember(wind_op, allowed_wind_ops)
                E.badinput('WIND_OP must be one of %s', strjoin(allowed_wind_ops, ', '));
            end
            
            if isempty(wind_speed)
                wind_speed = ask_number('Enter the wind speed criterion in m/s', 'testfxn', @(x) isscalar(x) && x >= 0, 'testmsg', 'Enter a scalar positive number');
            elseif ~isscalar(wind_speed) || ~isnumeric(wind_speed)
                E.badinput('WIND_SPEED must be a scalar number >= 0');
            end
        end
        
    end
    
end
