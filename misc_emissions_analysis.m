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
        
        function value = wrf_vcd_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.emis_wrf_dir, 'WRF-VCDs');
        end
        
        function filename = avg_file_name(year_in, days_of_week)
            years_str = strjoin(sprintfmulti('%d', year_in),'_');
            filename = sprintf('Summer_avg_%s_%s.mat', years_str, days_of_week);
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
        
        function locs = read_locs_file()
            locs = read_loc_spreadsheet();
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
        
        function [start_dates, end_dates, time_period] = select_start_end_dates(time_period)
            E = JLLErrors;
            % Eventually these should change to cell arrays once the line
            % density functions are set up to take non-contiguous time
            % periods
            if nargin < 1 || isempty(time_period)
                time_period = ask_multichoice('Which time period to use?', {'beginning (2007-09)', 'end (2012-14)'}, 'list', true);
                time_period = strsplit(time_period, ' ');
                time_period = time_period{1};
            end
            
            start_month = 4;
            end_month = 9;
            
            if strcmpi(time_period, 'beginning')
                start_dates = {datenum(2007, start_month, 1), datenum(2008, start_month, 1), datenum(2009, start_month, 1)};
                end_dates = {eomdate(2007, end_month), eomdate(2008, end_month), eomdate(2009, end_month)};
            elseif strcmpi(time_period, 'end')
                start_dates = {datenum(2012, start_month, 1), datenum(2013, start_month, 1), datenum(2014, start_month, 1)};
                end_dates = {eomdate(2012, end_month), eomdate(2013, end_month), eomdate(2014, end_month)};
            else
                E.badinput('TIME_PERIOD "%s" not recognized', time_period);
            end
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
        
        function [xx, radius] = find_indices_in_radius_around_loc(loc, lon, lat)
            % The Radius field is a carry over from Russell et al. 2012.
            % Now we use the boxes defined for the EMG fitting for
            % consistency.
            radius = mean(loc.BoxSize(3:4));
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
            
            % As of 13 Feb 2018, all of the fits files are for the first 70
            % locations, which misses Valmy, but I can fix that when I
            % rerun for the 3 year data anyway. Thus we only need to
            % override the default value if the fits files don't cover
            % sites 1-70.
            first_locs = load(misc_emissions_analysis.fits_file_name(first_dates_st(1), first_dates_end(end), use_wrf, file_loc_inds, first_weekdays, fit_type));
            second_locs = load(misc_emissions_analysis.fits_file_name(second_dates_st(1), second_dates_end(end), use_wrf, file_loc_inds, second_weekdays, fit_type));
            
            first_locs.locs = misc_emissions_analysis.cutdown_locs_by_index(first_locs.locs, user_loc_inds);
            second_locs.locs = misc_emissions_analysis.cutdown_locs_by_index(second_locs.locs, user_loc_inds);
            loc_names = {first_locs.locs.Location};
            loc_coords.lon = [first_locs.locs.Longitude]';
            loc_coords.lat = [first_locs.locs.Latitude]';
            
            % Collect the emissions and lifetimes differences into two
            % n-by-2 arrays. Also get some metrics of the goodness of fit
            % and total mass
            
            additional_fns = {'r2','a'};
            emis_tau_fns = fieldnames(first_locs.locs(1).emis_tau);
            all_fns = veccat(emis_tau_fns, additional_fns, 'column');
            default_mat = nan(numel(first_locs.locs), 2);
            for i_fn = 1:numel(all_fns)
                changes.(all_fns{i_fn}) = default_mat;
                for i_loc = 1:numel(first_locs.locs)
                    
                    first_value = find_substruct_field(first_locs.locs(i_loc), all_fns{i_fn});
                    second_value = find_substruct_field(second_locs.locs(i_loc), all_fns{i_fn});
                    if isempty(first_value)
                        first_value = nan;
                    end
                    if isempty(second_value)
                        second_value = nan;
                    end
                    
                    changes.(all_fns{i_fn})(i_loc, :) = [first_value, second_value];
                end
            end
            
            if include_vcds
                changes.vcds = default_mat;
                changes.vcds(:,1) = misc_emissions_analysis.avg_vcds_around_loc(first_locs.locs, first_time_period, first_weekdays);
                changes.vcds(:,2) = misc_emissions_analysis.avg_vcds_around_loc(second_locs.locs, second_time_period, second_weekdays);
            end
            
            changes.Location = {first_locs.locs.Location}';
            changes.ShortName = {first_locs.locs.ShortName}';
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
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generation methods %
        %%%%%%%%%%%%%%%%%%%%%%
        function make_summer_averages(varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('avg_year',[]);
            p.addParameter('days_of_week','');
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_year = pout.avg_year;
            days_of_week = pout.days_of_week;
            
            if isempty(avg_year)
                avg_year = ask_number('Enter the year (or years separated by a space) to do a summer average for', 'testfxn', @(x) all(x >= 2005 & x <= 2015), 'testmsg', 'Year(s) must be between 2005 and 2015');
            elseif ~isnumeric(avg_year) || any(avg_year < 2005 | avg_year > 2015)
                E.badinput('AVG_YEAR must be a numeric vector with values between 2005 and 2015');
            end
            
            days_of_week = misc_emissions_analysis.choose_days_of_week(days_of_week);
            
            start_date = cell(size(avg_year));
            end_date = cell(size(avg_year));
            for a=1:numel(avg_year)
                start_date{a} = datenum(avg_year(a), 4, 1);
                end_date{a} = datenum(avg_year(a), 9, 30);
            end
            
            % Make the monthly profile product average, then try to make
            % the daily one. If there's no data, it will return a NaN
            common_opts = {'DEBUG_LEVEL', 1, 'dayofweek', days_of_week};
            [monthly.no2, monthly.lon, monthly.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'monthly', common_opts{:});
            [daily.no2, daily.lon, daily.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'daily', common_opts{:});
            
            save_name = misc_emissions_analysis.avg_file_name(avg_year, days_of_week);
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
            weight_wind_dirs = pout.weight_wind_dirs;
            
            if ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('The parameter "loc_indicies" must be a numeric array with all values >= 1')
            end
            
            if (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('The parameter "do_overwrite" must be a scalar logical or number')
            end
            
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
 
            % Find the list of BEHR files between the start and end dates
            [behr_files, behr_dir] = list_behr_files(start_date, end_date,'daily','all');
            wind_reject_field = '';
            if wrf_bool
                for i_file = 1:numel(behr_files)
                    behr_files(i_file).name = strrep(behr_files(i_file).name, 'OMI', 'WRF');
                end
                behr_dir = misc_emissions_analysis.wrf_vcd_dir;
                wind_reject_field = 'WRFWindRejects';
            end
            % If we're doing line densities by sector, then we don't want
            % to reject any wind directions. We also don't want to reject
            % wind directions if explicitly told not to.
            filter_by_wind_dir = use_wind_rejects && ~by_sectors;
            if ~filter_by_wind_dir
                wind_reject_field = 'none';
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
                E.callError('date_mismatch', 'Dates in winds file (%s) do not match required (%s to %s, %d dates)', winds_file, datestr(start_date(1)), datestr(end_date(end)), numel(check_dvec));
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
            
            parfor a=1:numel(winds.locs)
                opt_args = {};
                
                box_size = winds_locs_distributed(a).BoxSize;
                if any(isnan(box_size))
                    warning('NaN detected in box size. Setting to default %s for location %s', mat2str(default_box), winds_locs_distributed(a).ShortName)
                    box_size = default_box;
                end
                
                if weight_wind_dirs && ~by_sectors
                    opt_args = veccat(opt_args, {'wind_dir_weights', wind_dir_weights{a}, 'wind_weights_bins', wind_dir_edges{a}});
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
        
        function make_avg_wrf_lifetimes(varargin)
            % This will compute the average WRF lifetimes vs. HNO3, ANs,
            % and total for WRF files across the specified time periods.
            % Since RO2 species aren't stored in the output, we need to
            % calculate them assuming steady state RO2 concentrations,
            % i.e. P(RO2) = L(RO2).
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
            
            p.parse(varargin{:});
            pout = p.Results;
            loc_inds = pout.loc_indicies;
            
            [beg_start_date, beg_end_date] = misc_emissions_analysis.select_start_end_dates('beginning');
            [end_start_date, end_end_date] = misc_emissions_analysis.select_start_end_dates('end');
            
            fit_type = misc_emissions_analysis.get_fit_type_interactive();
            
            % Load the weekday, weekend, and all day fits files. Cut out
            % the memory intensive parts of the line density sub structure
            % to save memory.
            [use_wrf, loc_inds_wrf] = misc_emissions_analysis.ask_to_use_wrf();
            if isempty(loc_inds)
                loc_inds = loc_inds_wrf;
            end
            
            days_of_week = {'UMTWRFS', 'TWRF', 'US'};
            loc_prototype = struct('Location', '', 'x', [], 'linedens', [], 'emgfit', [], 'r2', []);
            beg_locs = repmat(loc_prototype, numel(loc_inds), numel(days_of_week));
            end_locs = repmat(loc_prototype, numel(loc_inds), numel(days_of_week));
            
            for i_dow = 1:numel(days_of_week)
                beg_fits = load(misc_emissions_analysis.fits_file_name(beg_start_date, beg_end_date, use_wrf, loc_inds, days_of_week{i_dow}, fit_type));
                end_fits = load(misc_emissions_analysis.fits_file_name(end_start_date, end_end_date, use_wrf, loc_inds, days_of_week{i_dow}, fit_type));
                for i_loc = 1:numel(beg_fits.locs)
                    beg_locs(i_loc, i_dow) = copy_structure_fields(beg_fits.locs(i_loc), loc_prototype, 'substructs');
                    end_locs(i_loc, i_dow) = copy_structure_fields(end_fits.locs(i_loc), loc_prototype, 'substructs');
                end
            end
            
            % Now we can actually do the plotting. Make a list of the
            % available locations, then ask which one to plot until we quit
            locs_list = {beg_locs(:,1).Location};
            colors = {[0.5 0.5 0.5], 'c', [1 0.5 0];...
                      'k'          , 'b', 'r'};
            while true
                loc_ind = ask_multichoice('Plot which location?', locs_list, 'list', true, 'index', true);
                beg_fig = figure;
                l = gobjects(size(beg_locs,2),1);
                for i_wkday = 1:size(beg_locs,2)
                    l(i_wkday) = line(beg_locs(loc_ind, i_wkday).x, beg_locs(loc_ind, i_wkday).linedens, 'color', colors{1, i_wkday}, 'marker', 'o', 'linestyle', 'none');
                    if ~isempty(beg_locs(loc_ind, i_wkday).emgfit)
                        line(beg_locs(loc_ind, i_wkday).x, beg_locs(loc_ind, i_wkday).emgfit, 'color', colors{2, i_wkday}, 'linestyle', '--');
                    end
                end
                r2_array = {beg_locs(loc_ind, :).r2};
                xx = iscellcontents(r2_array, @isempty);
                r2_array(xx) = {nan};
                legend(l, sprintfmulti('%s (R^2 = %.2f)', days_of_week, r2_array));
                title(sprintf('Beginning, %s', beg_locs(loc_ind, 1).Location)); 
                
                end_fig = figure;
                l = gobjects(size(beg_locs,2),1);
                for i_wkday = 1:size(beg_locs,2)
                    l(i_wkday) = line(end_locs(loc_ind, i_wkday).x, end_locs(loc_ind, i_wkday).linedens, 'color', colors{1, i_wkday}, 'marker', 'o', 'linestyle', 'none');
                    if ~isempty(end_locs(loc_ind, i_wkday).emgfit)
                        line(end_locs(loc_ind, i_wkday).x, end_locs(loc_ind, i_wkday).emgfit, 'color', colors{2, i_wkday}, 'linestyle', '--');
                    end
                end
                r2_array = {end_locs(loc_ind, :).r2};
                xx = iscellcontents(r2_array, @isempty);
                r2_array(xx) = {nan};
                legend(l, sprintfmulti('%s (R^2 = %.2f)', days_of_week, r2_array));
                title(sprintf('End, %s', end_locs(loc_ind, 1).Location));
                tilefigs;
                
                input('Press ENTER to continue','s');
                
                close([beg_fig, end_fig]);
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
        
        function plot_lifetime_vs_mass(varargin)
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
            p.addParameter('use_wrf', nan);
            p.addParameter('single_plot', nan);
            p.addParameter('days_of_week', '');
            
            p.parse(varargin{:});
            pout = p.Results;
                        
            loc_inds = pout.loc_inds;
            if isnan(loc_inds)
                [loc_inds, file_loc_inds] = misc_emissions_analysis.get_loc_inds_interactive();
            end
            
            mass_value = pout.mass_value;
            allowed_mass_vals = {'a','vcds'};
            if isempty(mass_value)
                mass_value = ask_multichoice('Which quantity to use for mass of NOx?', allowed_mass_vals, 'list', true);
            elseif ~ismember(mass_value, allowed_mass_vals)
                E.badinput('MASS_VALUE must be one of: %s', strjoin(allowed_mass_vals, ', '));
            end
            
            use_wrf = pout.use_wrf;
            [use_wrf, loc_inds, ~, file_loc_inds] = misc_emissions_analysis.ask_to_use_wrf(use_wrf, loc_inds, file_loc_inds);
            
            single_plot_bool = pout.single_plot;
            if isnan(single_plot_bool)
                single_plot_bool = ask_yn('Plot all locations on a single plot?');
            elseif ~isscalar(single_plot_bool) || ~islogical(single_plot_bool)
                E.badinput('"single_plot" must be a scalar logical');
            end
            
            days_of_week = pout.days_of_week; % only used in single plot mode
            if single_plot_bool
                allowed_dows = {'UMTWRFS','TWRF','US'};
                if isempty(days_of_week)
                    days_of_week = ask_multiselect('Choose which day-of-week subsets to include', allowed_dows);
                elseif ~ismember(days_of_week, allowed_dows) && ~strcmpi(days_of_week, 'all')
                    E.badinput('"days_of_week must be one of: "%s", or "all"', strjoin(allowed_dows, '", "'));
                end
                if strcmpi(days_of_week, 'all')
                    days_of_week = allowed_dows;
                elseif ischar(days_of_week)
                    days_of_week = {days_of_week};
                end
            end
            
            fit_type = misc_emissions_analysis.get_fit_type_interactive(pout.fit_type);
            
            locs = misc_emissions_analysis.read_locs_file();
            
            vcds_bool = strcmpi(mass_value, 'vcds');
            decadal_changes = misc_emissions_analysis.collect_changes('beginning', 'end', 'UMTWRFS', 'UMTWRFS', 'loc_inds', loc_inds, 'file_loc_inds', file_loc_inds, 'use_wrf', use_wrf, 'include_vcds', vcds_bool, 'fit_type', fit_type);
            beginning_changes = misc_emissions_analysis.collect_changes('beginning', 'beginning', 'TWRF', 'US', 'loc_inds', loc_inds, 'file_loc_inds', file_loc_inds, 'use_wrf', use_wrf, 'include_vcds', vcds_bool, 'fit_type', fit_type);
            end_changes = misc_emissions_analysis.collect_changes('end', 'end', 'TWRF', 'US', 'loc_inds', loc_inds, 'file_loc_inds', file_loc_inds, 'use_wrf', use_wrf, 'include_vcds', vcds_bool, 'fit_type', fit_type);
            
            
            decadal_style = struct('marker', {'o','x'}, 'linestyle', 'none', 'color', {'b','r'},'markersize',10);
            beginning_style = struct('marker', {'^','*'}, 'linestyle', 'none', 'color', {'c','m'},'markersize',10);
            end_style = struct('marker', {'^','*'}, 'linestyle', 'none', 'color', {[1 0.5 0], [0.5 0 0.5]},'markersize',10);
            
            if vcds_bool
                x_label_str = 'Avg. NO_2 VCD (molec. cm^2)';
            else
                x_label_str = 'a (mol NO_2)';
            end
            
            if ~single_plot_bool
                for i_loc = 1:numel(loc_inds)
                    l = gobjects(6,1);
                    figure;
                    ax = gca;
                    l(1:2) = plot_changes(decadal_changes.(mass_value)(i_loc,:), decadal_changes.tau(i_loc,:), 'group_fmts', decadal_style, 'parent', ax);
                    l(3:4) = plot_changes(beginning_changes.(mass_value)(i_loc,:), beginning_changes.tau(i_loc,:), 'group_fmts', beginning_style, 'parent', ax);
                    l(5:6) = plot_changes(end_changes.(mass_value)(i_loc,:), end_changes.tau(i_loc,:), 'group_fmts', end_style, 'parent', ax);
                    legend(l, {'2005-07 UMTWRFS','2012-13 UMTWRFS','2005-07 TWRF','2005-07 SU','2012-13 TWRF','2012-13 SU'});
                    set(ax,'fontsize',16);
                    xlabel(x_label_str);
                    ylabel('\tau (hours)');
                    title(locs(loc_inds(i_loc)).Location);
                end
            else
                l = gobjects(0);
                legend_cell = {};
                figure;
                if ismember('UMTWRFS', days_of_week)
                    l(end+1) = line(decadal_changes.(mass_value)(:,1), decadal_changes.tau(:,1), decadal_style(1));
                    l(end+1) = line(decadal_changes.(mass_value)(:,2), decadal_changes.tau(:,2), decadal_style(2));
                    legend_cell = veccat(legend_cell, {'UMTWRFS 2005-2007', 'UMTWRFS 2012-2013'});
                end
                if ismember('TWRF', days_of_week)
                    l(end+1) = line(beginning_changes.(mass_value)(:,1), beginning_changes.tau(:,1), beginning_style(1));
                    l(end+1) = line(end_changes.(mass_value)(:,1), end_changes.tau(:,1), end_style(1));
                    legend_cell = veccat(legend_cell, {'TWRF 2005-2007', 'TWRF 2012-2013'});
                end
                if ismember('US', days_of_week)
                    l(end+1) = line(beginning_changes.(mass_value)(:,2), beginning_changes.tau(:,2), beginning_style(2));
                    l(end+1) = line(end_changes.(mass_value)(:,2), end_changes.tau(:,2), end_style(2));
                    legend_cell = veccat(legend_cell, {'US 2005-2007', 'US 2012-2013'});
                end
                
                legend(l(:), legend_cell);
                set(gca,'xscale','log','fontsize',14);
                xlabel(x_label_str);
                ylabel('\tau (hours)');
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
                    E.notimplemented('Converting winds from vector-type definition to met definition (from N is 0 deg)')
                else
                    north_angle = 90;
                    east_angle = 0;
                    fprintf('!!! NOTE: wind rose plotting directions winds are blowing TOWARDS !!!\n');
                end
                
                % Wind rose automatically creates a new figure
                WindRose(wind_dir, wind_vel, 'anglenorth', north_angle, 'angleeast', east_angle);
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
        
        function avg_vcds = avg_vcds_around_loc(locs, time_period, days_of_week)
            [start_dates, end_dates] = misc_emissions_analysis.select_start_end_dates(time_period);
            time_period_years = unique(cellfun(@year, veccat(start_dates, end_dates)));
            VCDs = load(misc_emissions_analysis.avg_file_name(time_period_years, days_of_week));
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
                xx_radius = misc_emissions_analysis.find_indices_in_radius_around_loc(locs(i_loc), VCDs.daily.lon, VCDs.daily.lat);
                avg_vcds(i_loc) = nanmean(VCDs.daily.no2(xx_radius));
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
            
            difference_is_significant = false(size(values,1),1);
            
            % Calculate whether the difference is significant at the 95%
            % confidence level.
            for i_chng = 1:size(values,1)
                % For SSR (sum of squared residuals), since SD = sqrt( SSR
                % / n_dofs ) calculate SSR from number of degrees of
                % freedom and square root. For the numbers of measurements,
                % we need to add back in the 5 degrees of freedom removed
                % in accounting for the fitting parameters.
                %
                % 23 Apr 2018 - strictly we could calculate the sum of
                % squared residuals by taking sum((emgfit - linedens).^2),
                % whereas doing sd^2 * n_dof includes the variance of the
                % fit in there. However, I'm not sure which is better for
                % the t-test.
                error('Did you figure out if the SSR should be sum(emgfit - linedens) or sd^2 * n?');
                [~, ~, difference_is_significant(i_chng)] = two_sample_t_test(values(i_chng, 1), value_sds(i_chng, 1).^2 .* value_dofs(i_chng, 1), value_dofs(i_chng, 1) + 5,...
                    values(i_chng, 2), value_sds(i_chng, 2).^2 .* value_dofs(i_chng, 2), value_dofs(i_chng, 2) + 5,...
                    sum(value_dofs(i_chng,:)), 'confidence', 0.95);
            end
            
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
        
        function locs = cutdown_locs_by_index(locs, loc_inds)
            % Cuts down the structure LOCS to the locations specified by
            % LOCS_INDS by matching the site names from LOCS against the
            % location names in the structure read in from the spreadsheet.
            % If LOCS_INDS is an empty array, do not cut the locations down
            % at all.
            
            E = JLLErrors;
            
            if isempty(loc_inds)
                return 
            end
            
            locs_ss = misc_emissions_analysis.read_locs_file();
            ss_names = {locs_ss(loc_inds).Location};
            in_names = {locs.Location};
            xx = ismember(in_names, ss_names);
            if sum(xx) ~= numel(loc_inds)
                E.callError('loc_lookup_error', 'A different number of locations was found in LOCS that specified by LOC_INDS');
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
            % 26 Feb 2018: WRF line densities only completed for Chicago,
            % so if using WRF line densities, we need to specify location
            % indicies = 9 for the file name.
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
                loc_inds = 9; % currently WRF data only available for Chicago
                wind_rej_field = misc_emissions_analysis.wind_reject_field_wrf;
                file_loc_inds = loc_inds;
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
