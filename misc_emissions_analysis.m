classdef misc_emissions_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        % These are the wind speed division and whether to use fast or slow
        % wind speeds to generate the rotated line densities.
        em_wind_spd = 3;
        em_wind_mode = 'fast';
        
    end
    
    properties(Constant = true, Access = protected)
        % This is used to help check if all the required Git repos are in
        % the proper state. The first time any method calls the Git
        % verification method, if it passes, this is set to true so that it
        % doesn't need to be checked again.
        git_check_complete = false;
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
        
        function value = emis_wrf_dir
            value = misc_emissions_analysis.subdir_prep(misc_emissions_analysis.workspace_dir, 'WRFData');
        end
        
        function filename = avg_file_name(year_in)
            years_str = strjoin(sprintfmulti('%d', year_in),'_');
            filename = sprintf('Summer_avg_%s.mat', years_str);
            filename = fullfile(misc_emissions_analysis.avg_save_dir, filename);
        end
        
        function filename = winds_file_name(start_date, end_date)
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('site_winds_%sto%s.mat', datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(end), 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.site_info_dir, filename);
        end
        
        function filename = line_density_file_name(start_date, end_date, by_sectors, all_winds_bool)
            if by_sectors
                sectors_string = 'sectors';
            else
                sectors_string = 'rotated';
            end
            if all_winds_bool
                winds_string = '_allwinds';
            else
                winds_string = '';
            end
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('site_%s%s_no2_%sto%s.mat', sectors_string, winds_string, datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(2), 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.line_density_dir, filename);
        end
        
        function filename = fits_file_name(start_date, end_date, is_subset)
            if is_subset
                subset_str = '_subset';
            else
                subset_str = '';
            end
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            filename = sprintf('site_emg_fits%s_%sto%s.mat', subset_str, datestr(start_date(1), 'yyyy-mm-dd'), datestr(end_date(2), 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.line_density_dir, filename);
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
                % to handle discontinuous date ranges.
                G = GitChecker;
                G.addReqCommits(behr_paths.behr_utils, 'd524710e');
                G.addReqCommits(behr_paths.utils, '9af1577');
                % checkState() by default will error if the repositories
                % are not in the correct state.
                G.checkState()
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
        
        function [start_dates, end_dates] = select_start_end_dates(time_period)
            E = JLLErrors;
            % Eventually these should change to cell arrays once the line
            % density functions are set up to take non-contiguous time
            % periods
            if nargin < 1 || isempty(time_period)
                time_period = ask_multichoice('Which time period to use?', {'beginning (2005, 2007)', 'end (2012-13)'}, 'list', true);
                time_period = strsplit(time_period, ' ');
                time_period = time_period{1};
            end
            
            start_month = 4;
            end_month = 9;
            
            if strcmpi(time_period, 'beginning')
                start_dates = {datenum(2005, start_month, 1), datenum(2007, start_month, 1)};
                end_dates = {eomdate(2005, end_month), eomdate(2007, end_month)};
            elseif strcmpi(time_period, 'end')
                start_dates = {datenum(2012, start_month, 1), datenum(2013, start_month, 1)};
                end_dates = {eomdate(2013, end_month), eomdate(2013, end_month)};
            else
                E.badinput('TIME_PERIOD "%s" not recognized', time_period);
            end
        end
        
        function [xx, yy] = find_loc_indices(loc, lon, lat, radius)
            % LOC must be a scalar element of the locations structure, LON
            % and LAT must be 2D arrays of longitude and latitude
            % coordinates for an NO2 average or similar 2D field. RADIUS
            % must be a scalar number of grid cells in each direction to
            % get. If omitted, defaults to 0.
            E = JLLErrors;
            
            if ~exist('radius', 'var')
                radius = 0;
            end
            
            r = sqrt((lon - loc.Longitude).^2 + (lat - loc.Latitude).^2);
            [~, i_min] = min(r(:));
            [xx, yy] = ind2sub(size(lon), i_min); 
            
            xx = (xx - radius):(xx + radius);
            yy = (yy - radius):(yy + radius);
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
        
        function wind_logical = set_wind_conditions(location, speed_cutoff, fast_or_slow)
            E = JLLErrors;
            if ~isstruct(location) || ~isscalar(location) || any(~isfield(location, {'ShortName', 'WindDir', 'WindSpeed'}))
                E.badinput('LOCATION must be a scalar structure with fields "ShortName", "WindDir", and "WindSpeed"')
            end
            
            if ~isnumeric(speed_cutoff) || ~isscalar(speed_cutoff) || speed_cutoff < 0
                E.badinput('SPEED_CUTOFF must be a scalar, positive number')
            end
            
            allowed_fast_slow = {'fast', 'slow', 'none'};
            if ~ismember(fast_or_slow, allowed_fast_slow)
                E.badinput('FAST_OR_SLOW must be one of: %s', strjoin(allowed_fast_slow));
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
            for a=1:size(location.WindRejects,1)
                wind_dir_range = location.WindRejects(a,:);
                if wind_dir_range(1) < wind_dir_range(2)
                    xx = location.WindDir >= wind_dir_range(1) & location.WindDir < wind_dir_range(2);
                else
                    xx = location.WindDir >= wind_dir_range(1) | location.WindDir < wind_dir_range(2);
                end
                
                wind_logical(xx) = false;
            end
            
            % Handle wind speed filtering here
            if strcmpi(fast_or_slow, 'fast')
                wind_logical(location.WindSpeed < speed_cutoff) = false;
            elseif strcmpi(fast_or_slow, 'slow')
                wind_logical(location.WindSpeed >= speed_cutoff) = false;
            elseif ~strcmpi(fast_or_slow, 'none')
                E.badinput('FAST_OR_SLOW must be one of: "fast", "slow", or "none"');
            end
        end
        
        function emis_tau = calculate_emission_lifetime(line_dens_struct, fit_struct, wind_speed_vector)
            % First we need to compute the total uncertainty in the
            % parameters that accounts for uncertainty in the NO2 VCDs,
            % across wind integration distance, choice of wind fields, etc.
            param_uncert = calc_fit_param_uncert(fit_struct.ffit, fit_struct.param_stats.percent_ci95/100, line_dens_struct.num_valid_obs);
            
            % Then use these uncertainties to calculate the emissions and
            % lifetime and their uncertainties
            [emis_tau.emis, emis_tau.emis_uncert, emis_tau.tau, emis_tau.tau_uncert] = compute_emg_emis_tau(fit_struct.ffit.a, param_uncert(1), fit_struct.ffit.x_0, param_uncert(2), 'vec', wind_speed_vector, 'emissions_type', 'no');
        end
        
        function [nei_no, nei_lon, nei_lat] = load_nei_by_year(nei_year)
            % Make the input path
            tmp_path = find_wrf_path('us','daily',datenum(nei_year,1,1));
            path_parts = strsplit(tmp_path, '/');
            
            % The path on my computer is something like '/Volumes/share-wrfN/...' 
            % and we just want to get which network drive it should be on.
            % The first three parts of the split path should be an empty
            % string, Volumes, and the share. This will put a / at the
            % beginning
            wrf_share = strjoin(path_parts(1:3), '/');
            
            % Whichever share it's on, it should be in a consistent path
            % there - except for 2012. I was having trouble getting the
            % full year's inputs to prepare, so I had to split it into two
            % 6 month periods. The NEI emissions are the same in both, so
            % we can just pick one.
            inputs_path = fullfile(wrf_share, 'Inputs', num2str(nei_year), 'IC-BC-Emis');
            if nei_year == 2012
                inputs_path = fullfile(inputs_path, 'Months01-06');
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
            % Convert from mol NO / km^2 / hr to Mg NO / hr: molar mass of NO =
            % 30.006 g / mol = 30.006e-6 Mg / mol
            nei_no = nei_no .* grid_area .* 30.06e-6; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generation methods %
        %%%%%%%%%%%%%%%%%%%%%%
        function make_summer_averages(avg_year)
            if ~exist('avg_year', 'var')
                avg_year = ask_number('Enter the year (or years separated by a space) to do a summer average for', 'testfxn', @(x) all(x >= 2005 & x <= 2015), 'testmsg', 'Year(s) must be between 2005 and 2015');
            end
            
            start_date = cell(size(avg_year));
            end_date = cell(size(avg_year));
            for a=1:numel(avg_year)
                start_date{a} = datenum(avg_year(a), 4, 1);
                end_date{a} = datenum(avg_year(a), 9, 30);
            end
            
            % Make the monthly profile product average, then try to make
            % the daily one. If there's no data, it will return a NaN
            common_opts = {'DEBUG_LEVEL', 1};
            [monthly.no2, monthly.lon, monthly.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'monthly', common_opts{:});
            [daily.no2, daily.lon, daily.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'daily', common_opts{:});
            
            save_name = misc_emissions_analysis.avg_file_name(avg_year);
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
                if month(dvec(d)) == last_month && year(dvec(d)) == last_year
                    fprintf('     Using existing list of WRF files\n');
                    wrf_files = misc_emissions_analysis.closest_wrf_file_in_time(dvec(d), all_months_wrf_files);
                else
                    fprintf('     New month: need to get the directory listing\n');
                    [wrf_files, all_months_wrf_files] = misc_emissions_analysis.closest_wrf_file_in_time(dvec(d));
                    last_month = month(dvec(d));
                    last_year = year(dvec(d));
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
                        [xx,yy] = misc_emissions_analysis.find_loc_indices(locs(l), wrf_lon, wrf_lat, 1);
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
            
            write_date = datestr(now); %#ok<NASGU>
            
            save(save_file, 'locs', 'dvec', 'write_date');
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
            misc_emissions_analysis.make_line_densities(true, varargin{:});
        end
        
        function make_line_densities(by_sectors, varargin)
            
            if ~islogical(by_sectors) || ~isscalar(by_sectors)
                E.badinput('BY_SECTORS must be a scalar logical')
            end
            
            p = advInputParser;
            p.addParameter('time_period', '');
            p.addParameter('loc_indices', []);
            p.addParameter('do_overwrite', -1);
            p.addFlag('all_winds');
            
            p.parse(varargin{:});
            pout = p.AdvResults;
            
            misc_emissions_analysis.verify_git_state();
            
            time_period = pout.time_period;
            loc_indicies = pout.loc_indices;
            do_overwrite = pout.do_overwrite;
            all_winds = pout.all_winds;
            
            if ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('The parameter "loc_indicies" must be a numeric array with all values >= 1')
            end
            
            if (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('The parameter "do_overwrite" must be a scalar logical or number')
            end
            
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
            
            % If overwrite not given and the save file exists, ask to
            % overwrite. Otherwise, only overwrite if instructed.
            save_name = misc_emissions_analysis.line_density_file_name(start_date, end_date, by_sectors, all_winds);
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
 
            % Find the list of BEHR files between the start and end dates
            [behr_files, behr_dir] = list_behr_files(start_date, end_date,'daily');
            
            % Load the winds file with up-to-date box sizes and wind
            % sectors to reject
            winds = misc_emissions_analysis.load_winds_file(start_date, end_date);
            
            % Check that the dates match up with what we're expecting (it
            % should because we load the file with those dates)
            check_dvec = misc_emissions_analysis.make_datevec(start_date, end_date);
            if ~isequal(check_dvec, winds.dvec)
                E.callError('date_mismatch', 'Dates in winds file (%s) do not match required (%s to %s, %d dates)', winds_file, datestr(start_date(1)), datestr(end_date(end)), numel(check_dvec));
            end
            
            if ~isempty(loc_indicies)
                winds.locs = winds.locs(loc_indicies);
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
            
            parfor a=1:numel(winds.locs)
                fprintf('Calculating sector line densities for %s\n', winds_locs_distributed(a).ShortName);
                
                box_size = winds_locs_distributed(a).BoxSize;
                if any(isnan(box_size))
                    warning('NaN detected in box size. Setting to default %s for location %s', mat2str(default_box), winds_locs_distributed(a).ShortName)
                    box_size = default_box;
                end
                
                % Choose slow wind days for this - our goal is to find out
                % which directions have downwind sources that will confound the
                % fitting
                if by_sectors
                    % Place holder - I only use the sectors code to find
                    % directions that are good for the line density
                    % analysis, so I look at slow winds to find directions
                    % that have downwind sources (like a manual version of
                    % Liu et al. 2016).
                    %
                    % Alternately, since I really just want the best maps I
                    % can get to look for downwind interferences, I can
                    % just use all wind speeds
                    if all_winds
                        wind_logical = true(size(winds_locs_distributed(a).WindSpeed));
                    else
                        wind_logical = winds_locs_distributed(a).WindSpeed < 3;
                    end
                    % Call the sector division algorithm
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density_sectors(behr_dir, behr_files, winds_locs_distributed(a).Longitude, winds_locs_distributed(a).Latitude, winds_locs_distributed(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size);
                else
                    wind_logical = misc_emissions_analysis.set_wind_conditions(winds_locs_distributed(a), misc_emissions_analysis.em_wind_spd, misc_emissions_analysis.em_wind_mode);
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density(behr_dir, behr_files, winds_locs_distributed(a).Longitude, winds_locs_distributed(a).Latitude, winds_locs_distributed(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size);
                end
                
                
                %winds.locs(a).no2_sectors = no2;
            end
            
            for a=1:numel(winds.locs)
                winds.locs(a).no2_sectors = no2(a);
            end
            
            locs = winds.locs;
            dvec = winds.dvec;
            write_date = datestr(now);
            
            save(save_name, '-v7.3', 'locs', 'dvec', 'write_date');
        end
        
        function locs = make_emg_fits(time_period, loc_indicies, add_nei, do_overwrite)
            E = JLLErrors;
            
            if ~exist('time_period', 'var')
                time_period = '';
            end

            if ~exist('loc_indicies', 'var')
                loc_indicies = [];
            elseif ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('LOC_INDICIES must be a numeric array with all values >= 1')
            end
            
            if ~exist('add_nei', 'var')
                add_nei = true;
            elseif ~isscalar(add_nei) || (~islogical(add_nei) && ~isnumeric(add_nei))
                E.badinput('ADD_NEI must be a scalar logical or numeric value');
            end
            
            if ~exist('do_overwrite', 'var')
                do_overwrite = -1;
            elseif (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('DO_OVERWRITE must be a scalar logical or number')
            end
            
            [start_date, end_date] = misc_emissions_analysis.select_start_end_dates(time_period);
            % If overwrite not given and the save file exists, ask to
            % overwrite. Otherwise, only overwrite if instructed.
            save_name = misc_emissions_analysis.fits_file_name(start_date, end_date, ~isempty(loc_indicies));
            if exist(save_name, 'file')
                if do_overwrite < 0
                    if ~ask_yn(sprintf('%s exists. Overwrite?', save_name))
                        return
                    end
                elseif ~do_overwrite
                    return
                end
            end
            
            % Load the file with the line densities
            ldens_file = misc_emissions_analysis.line_density_file_name(start_date, end_date, false);
            line_densities = load(ldens_file);
            
            % Check that the dates match up with what we're expecting (it
            % should because we load the file with those dates)
            check_dvec = misc_emissions_analysis.make_datevec(start_date, end_date);
            if ~isequal(check_dvec, line_densities.dvec)
                E.callError('date_mismatch', 'Dates in winds file (%s) do not match required (%s to %s)', ldens_file, datestr(start_date(1)), datestr(end_date(end)));
            end
            
            if ~isempty(loc_indicies)
                locs = line_densities.locs(loc_indicies);
            else
                locs = line_densities.locs;
            end
            
            % Load the NEI data. Will need to get lat/lon from
            % wrfinput_d01, b/c the wrfchemi files don't include lat-lon.
            if add_nei
                nei_year = unique(year(check_dvec));
                [nei_avg_no, nei_lon, nei_lat] = misc_emissions_analysis.load_nei_by_year(nei_year);
            end
            % Specify even the default options so that if fit_line_density
            % changes, we know exactly what options we wanted.
            common_opts = {'fmincon_output', 'none', 'emgtype', 'lu', 'fittype', 'ssresid', 'nattempts', 20};
            
            for a=1:numel(locs)
                fprintf('Fitting %s\n', locs(a).ShortName);
                while true
                    [fit.ffit, fit.emgfit, fit.param_stats, fit.f0, fit.history, fit.fitresults] = fit_line_density(locs(a).no2_sectors.x, locs(a).no2_sectors.linedens, common_opts{:});
                    % Try this a second time - if it gives a different
                    % answer, we should re-run, since that suggests we
                    % didn't find the minimum one time.
                    ffit = fit_line_density(locs(a).no2_sectors.x, locs(a).no2_sectors.linedens, common_opts{:});
                    
                    % Check that the two are the same to within 1%
                    diff_tolerance = 0.01;
                    rdel = reldiff(struct2array(ffit), struct2array(fit.ffit));
                    if all(abs(rdel) < diff_tolerance)
                        break
                    else
                        fprintf('Fit results differ by > %f%% (%s vs %s); retrying\n', diff_tolerance*100, struct2string(fit.ffit), struct2string(fit));
                    end
                end
                locs(a).fit_info = fit;
                
                % Add the emissions and lifetime. Use the 95% confidence
                % intervals as the uncertainty. We need to restrict the
                % winds to what should have been used to calculate the line
                % densities.
                wind_logical = misc_emissions_analysis.set_wind_conditions(locs(a), misc_emissions_analysis.em_wind_spd, misc_emissions_analysis.em_wind_mode);
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
            end
            
            dvec = line_densities.dvec;
            write_date = datestr(now);
            
            save(save_name, '-v7.3', 'locs', 'dvec', 'write_date');
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % Plotting methods %
        %%%%%%%%%%%%%%%%%%%%
        
        function plot_site_summer_avg(loc_to_plot, plot_year, monthly_or_daily, plot_axis)
            E = JLLErrors;
            
            locs = misc_emissions_analysis.read_locs_file();
            loc_names = {locs.ShortName};
            if ~exist('loc_to_plot', 'var') || isempty(loc_to_plot)
                loc_to_plot = ask_multichoice('Which location to plot?', loc_names, 'list', true);
            elseif ~any(strcmpi(loc_to_plot, loc_names))
                E.badinput('LOC_TO_PLOT must be one of the shortname in the trend locations file');
            end
            
            i_loc = strcmpi(loc_names, loc_to_plot);
            
            avg_files = dir(fullfile(misc_emissions_analysis.avg_save_dir, '*.mat'));
            avg_files = {avg_files.name};
            if ~exist('plot_year', 'var')
                file_to_plot = ask_multichoice('Plot which avg. file?', avg_files, 'list', true);
                file_to_plot = fullfile(misc_emissions_analysis.avg_save_dir, file_to_plot);
            else
                if ~isnumeric(plot_year) || ~isscalar(plot_year)
                    E.badinput('PLOT_YEAR must be a scalar number')
                end
                file_to_plot = misc_emissions_analysis.avg_file_name(plot_year);
                if ~exist(file_to_plot, 'file')
                    E.badinput('No average file for year %d', plot_year);
                end
            end
            
            allowed_mod_strings = {'both', 'monthly', 'daily'};
            if ~exist('monthly_or_daily', 'var')
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
            [yy, xx] = misc_emissions_analysis.find_loc_indices(locs(i_loc), avgs.monthly.lon, avgs.monthly.lat, n_cells);
            
            loc_longrid = avgs.monthly.lon(yy,xx);
            loc_latgrid = avgs.monthly.lat(yy,xx);
            if ismember(monthly_or_daily, {'both', 'monthly'})
                loc_no2grid = avgs.monthly.no2(yy,xx);
                figure;
                pcolor(loc_longrid, loc_latgrid, loc_no2grid);
                line(locs(i_loc).Longitude, locs(i_loc).Latitude, 'marker', 'p', 'color', 'k', 'linestyle', 'none');
                colorbar;
                caxis(calc_plot_limits(loc_no2grid(:), 1e15, 'max', [0 Inf]));
                title(sprintf('%s - Monthly profiles', locs(i_loc).ShortName));
            end
            % If daily profile avgs are available too, plot them as well
            if ~isscalar(avgs.daily.no2) && ismember(monthly_or_daily, {'both', 'daily'}) % if no data, the no2 grid will just be a scalar NaN
                loc_no2grid = avgs.daily.no2(yy,xx);
                if ~exist('plot_axis','var')
                    figure; 
                    pcolor(loc_longrid, loc_latgrid, loc_no2grid);
                else
                    pcolor(plot_axis, loc_longrid, loc_latgrid, loc_no2grid);
                end
                line(locs(i_loc).Longitude, locs(i_loc).Latitude, 'marker', 'p', 'color', 'k', 'linestyle', 'none');
                colorbar;
                caxis(calc_plot_limits(loc_no2grid(:), 1e15, 'max', [0 Inf]));
                title(sprintf('%s - Daily profiles', locs(i_loc).ShortName));
            end
        end
        
        function plot_site_sectors_linedens(locs_to_plot, locs_dvec)
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
            if numel(locs_year) > 1
                E.notimplemented('Line densities across multiple years')
            end
            
            % Map the subplot index to the proper direction
            direction_names = {'NW','N','NE','W','vcds','E','SW','S','SE'};
            for a=1:numel(locs_to_plot)
                figure;
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
        
        function plot_sector_no2avg_with_boxes(locs_to_plot)
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
                figure;
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
            % For now, I'm just going to assume that we want to plot all
            % the locations in the subset emissions file. Later, I'll add
            % the ability to choose a subset of locations.
            allowed_time_periods = {'beginning','end','both'};
            plot_time_period = ask_multichoice('Which time period(s) to plot?', allowed_time_periods, 'list', true);
            plot_beginning = any(strcmpi(plot_time_period, {'beginning','both'}));
            plot_end = any(strcmpi(plot_time_period, {'end','both'}));
            
            % Loading both time periods isn't very hard, so we'll load both
            % and then just plot the one we want.
            [beg_start_date, beg_end_date] = misc_emissions_analysis.select_start_end_dates('beginning');
            [end_start_date, end_end_date] = misc_emissions_analysis.select_start_end_dates('end');
            D = load(misc_emissions_analysis.fits_file_name(beg_start_date, beg_end_date, true));
            beg_locs = D.locs;
            D = load(misc_emissions_analysis.fits_file_name(end_start_date, end_end_date, true));
            end_locs = D.locs;
            
            % Extract the site name, satellite, and NEI derived emissions
            names = cell(size(beg_locs));
            beg_sat_emis = nan(size(beg_locs));
            beg_sat_errors = nan(size(beg_locs));
            beg_nei_emis = nan(size(beg_locs));
            beg_nei_errors = nan(size(beg_locs));
            end_sat_emis = nan(size(end_locs));
            end_sat_errors = nan(size(end_locs));
            end_nei_emis = nan(size(end_locs));
            end_nei_errors = nan(size(end_locs));
            for a=1:numel(beg_locs)
                names{a} = beg_locs(a).ShortName;
                beg_sat_emis(a) = beg_locs(a).emis_tau.emis;
                beg_sat_errors(a) = beg_locs(a).emis_tau.emis_uncert;
                beg_nei_emis(a) = beg_locs(a).emis_tau.nei_emis;
                end_sat_emis(a) = end_locs(a).emis_tau.emis;
                end_sat_errors(a) = end_locs(a).emis_tau.emis_uncert;
                end_nei_emis(a) = end_locs(a).emis_tau.nei_emis;
            end
            
            plot_data = [];
            legend_strings = {};
            if plot_beginning
                plot_data = cat(2, plot_data, beg_sat_emis, beg_nei_emis);
                legend_strings = cat(2, legend_strings, {'BEHR (2005)', 'NEI (2005)'});
            end
            if plot_end
                plot_data = cat(2, plot_data, end_sat_emis, end_nei_emis);
                legend_strings = cat(2, legend_strings, {'BEHR (2012)', 'NEI (2012)'});
            end
            
            figure;
            bar(plot_data);
            ylabel('NO Emissions (Mg h^{-1})');
            %bar_errors([sat_emis, nei_emis], [sat_errors, nei_errors]);
            set(gca,'xticklabel',names,'fontsize',16,'ygrid','on','xtickLabelRotation',-30);
            legend(legend_strings{:});
        end
    end
    
end

