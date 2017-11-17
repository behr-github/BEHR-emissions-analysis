classdef misc_emissions_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)
        em_start_date = '2012-04-01';
        em_end_date = '2012-09-30';
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
        
        function filename = avg_file_name(year_in)
            filename = sprintf('Summer_avg_%d.mat', year_in);
            filename = fullfile(misc_emissions_analysis.avg_save_dir, filename);
        end
        
        function filename = winds_file_name(start_date, end_date)
            filename = sprintf('site_winds_%sto%s.mat', datestr(start_date, 'yyyy-mm-dd'), datestr(end_date, 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.site_info_dir, filename);
        end
        
        function filename = sectors_file_name(start_date, end_date, by_sectors)
            if by_sectors
                sectors_string = 'sectors';
            else
                sectors_string = 'rotated';
            end
            filename = sprintf('site_%s_no2_%sto%s.mat', sectors_string, datestr(start_date, 'yyyy-mm-dd'), datestr(end_date, 'yyyy-mm-dd'));
            filename = fullfile(misc_emissions_analysis.line_density_dir, filename);
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
        
        function locs = read_locs_file()
            locs = read_loc_spreadsheet();
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
        
        function wind_logical = set_wind_conditions(location, speed_cutoff, fast_or_slow, error_if_undef)
            E = JLLErrors;
            if ~isstruct(location) || ~isscalar(location) || any(~isfield(location, {'ShortName', 'WindDir', 'WindSpeed'}))
                E.badinput('LOCATION must be a scalar structure with fields "ShortName", "WindDir", and "WindSpeed"')
            end
            
            if ~isnumeric(speed_cutoff) || ~isscalar(speed_cutoff) || speed_cutoff < 0
                E.badinput('SPEED_CUTOFF must be a scalar, positive number')
            end
            
            allowed_fast_slow = {'fast', 'slow'};
            if ~ismember(fast_or_slow, allowed_fast_slow)
                E.badinput('FAST_OR_SLOW must be one of: %s', strjoin(fast_or_slow));
            end
            
            if ~exist('error_if_undef','var')
                error_if_undef = false;
            elseif ~islogical(error_if_undef) || ~isscalar(error_if_undef)
                E.badinput('ERROR_IF_UNDEF must be a scalar logical')
            end
            
            wind_logical = true(size(location.WindDir));
            % These can be used as constants to filter particular
            % directions, or one can of could specify your own specific
            % directions as degrees.
            dir_N = [67.5 112.5];
            dir_NE = [22.5 67.5];
            dir_E = [-22.5 22.5];
            dir_SE = [-67.2 22.5];
            dir_S = [-112.5 -67.5];
            dir_SW = [-152.5 -112.5];
            dir_W = [152.5 -152.5];
            dir_NW = [112.5 152.5];
            
            % Define how you want to filter wind directions for locations
            % based on their short name here.
            switch location.ShortName
                case 'Bakersfield'
                    % Using just 2012, Bakersfield is really hard to tell.
                    % There always seems to be a downwind source, so for
                    % now just going to try using all directions.
                case 'Chicago'
                    % The directions that have weird behavior just seem to
                    % be due to missing data, averaging all wind directions
                    % should help with that.
                case 'Fresno'
                    % This may just be sampling issues, i.e. the sectors
                    % line densities had too few days to work from
                    wind_logical(location.WindDir > dir_SW(1) & location.WindDir < dir_SW(2)) = false;
                    wind_logical(location.WindDir > dir_N(1) & location.WindDir < dir_N(2)) = false;
                case 'Houston'
                    % North is a little weird, but could probably be used
                    % with < 100 km downwind
                    wind_logical(location.WindDir > dir_N(1) & location.WindDir < dir_N(2)) = false;
                case 'Washington DC'
                    % pass
                otherwise
                    msg = sprintf('No wind direction filtering for %s defined', location.ShortName);
                    if error_if_undef
                        E.notimplemented(msg);
                    else
                        fprintf('%s\n',msg);
                    end
            end
            
            % Handle wind speed filtering here
            if strcmpi(fast_or_slow, 'fast')
                wind_logical(location.WindSpeed < speed_cutoff) = false;
            else
                wind_logical(location.WindSpeed >= speed_cutoff) = false;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Generation methods %
        %%%%%%%%%%%%%%%%%%%%%%
        function make_summer_averages(avg_year)
            if ~exist('avg_year', 'var')
                avg_year = ask_number('Enter the year to do a summer average for', 'testfxn', @(x) isscalar(x) && x >= 2005 && x <= 2015, 'testmsg', 'Year must be between 2005 and 2015');
            end
            
            start_date = datenum(avg_year, 4, 1);
            end_date = datenum(avg_year, 9, 30);
            
            % Make the monthly profile product average, then try to make
            % the daily one. If there's no data, it will return a NaN
            common_opts = {'rejectmode', 'behr', 'DEBUG_LEVEL', 1};
            [monthly.no2, monthly.lon, monthly.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'monthly', common_opts{:});
            [daily.no2, daily.lon, daily.lat] = behr_time_average(start_date, end_date, 'prof_mode', 'daily', common_opts{:});
            
            save_name = misc_emissions_analysis.avg_file_name(avg_year);
            save(save_name, 'monthly', 'daily');
        end
        
        function make_location_winds_file(start_date, end_date, overwrite)
            % As in Laughner, Zare, and Cohen (2016, ACP) we will calculate
            % wind direction by averaging over the first 5 WRF layers in a
            % 3x3 grid centered on each location.
            
            if ~exist('overwrite', 'var')
                overwrite = -1;
            end
            
            locs = misc_emissions_analysis.read_locs_file();
            start_date = validate_date(start_date);
            end_date = validate_date(end_date);
            dvec = start_date:end_date;
            
            % Check that the save file exists
            save_file = misc_emissions_analysis.winds_file_name;
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
        
        function make_line_densities(by_sectors, loc_indicies, do_overwrite)
            if ~islogical(by_sectors) || ~isscalar(by_sectors)
                E.badinput('BY_SECTORS must be a scalar logical')
            end
            
            if ~exist('loc_indicies', 'var')
                loc_indicies = [];
            elseif ~isnumeric(loc_indicies) || any(loc_indicies(:) < 1)
                E.badinput('LOC_INDICIES must be a numeric array with all values >= 1')
            end
            
            if ~exist('do_overwrite', 'var')
                do_overwrite = -1;
            elseif (~isnumeric(do_overwrite) && ~islogical(do_overwrite)) || ~isscalar(do_overwrite)
                E.badinput('DO_OVERWRITE must be a scalar logical or number')
            end
            
            start_date = misc_emissions_analysis.em_start_date;
            end_date = misc_emissions_analysis.em_end_date;
            
            % If overwrite not given and the save file exists, ask to
            % overwrite. Otherwise, only overwrite if instructed.
            save_name = misc_emissions_analysis.sectors_file_name(start_date, end_date, by_sectors);
            if exist(save_name, 'file')
                if do_overwrite < 0
                    if ~ask_yn('%s exists. Overwrite?')
                        return
                    end
                elseif ~do_overwrite
                    return
                end
            end
            
            % Find the list of BEHR files between the start and end dates
            [behr_files, behr_dir] = list_behr_files(start_date, end_date,'daily');
            
            % Load the winds file
            winds_file = misc_emissions_analysis.winds_file_name(start_date, end_date);
            winds = load(winds_file);
            
            % Check that the dates match up with what we're expecting (it
            % should because we load the file with those dates)
            check_dvec = datenum(start_date):datenum(end_date);
            if ~isequal(check_dvec, winds.dvec)
                E.callError('date_mismatch', 'Dates in winds file (%s) do not match required (%s to %s)', winds_file, datestr(start_date), datestr(end_date));
            end
            
            if ~isempty(loc_indicies)
                winds.locs = winds.locs(loc_indicies);
            end
            % Rural sites aren't going to be that interesting
            xx = strcmpi('Cities', {winds.locs.SiteType}) | strcmpi('PowerPlants', {winds.locs.SiteType});
            winds.locs(~xx) = [];
            box_size = [1 2 1 1];
            
            for a=1:numel(winds.locs)
                fprintf('Calculating sector line densities for %s\n', winds.locs(a).ShortName);
                
                % Choose slow wind days for this - our goal is to find out
                % which directions have downwind sources that will confound the
                % fitting
                if by_sectors
                    % Place holder - I only use the sectors code to find
                    % directions that are good for the line density
                    % analysis, so I look at slow winds to find directions
                    % that have downwind sources (like a manual version of
                    % Liu et al. 2016)
                    wind_logical = winds.locs(a).WindSpeed < 3;
                    % Call the sector division algorithm
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density_sectors(behr_dir, behr_files, winds.locs(a).Longitude, winds.locs(a).Latitude, winds.locs(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size);
                else
                    wind_logical = misc_emissions_analysis.set_wind_conditions(winds.locs(a), 3, 'fast');
                    [no2(a).x, no2(a).linedens, no2(a).linedens_std, no2(a).lon, no2(a).lat, no2(a).no2_mean, no2(a).no2_std, no2(a).num_valid_obs, no2(a).nox, no2(a).debug_cell] ...
                        = calc_line_density(behr_dir, behr_files, winds.locs(a).Longitude, winds.locs(a).Latitude, winds.locs(a).WindDir, wind_logical, 'interp', false, 'rel_box_corners', box_size);
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
        
        %%%%%%%%%%%%%%%%%%%%
        % Plotting methods %
        %%%%%%%%%%%%%%%%%%%%
        
        function plot_site_summer_avg(loc_to_plot, plot_year, monthly_or_daily, plot_axis)
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
    end
    
end

