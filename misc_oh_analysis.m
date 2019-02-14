classdef misc_oh_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property-like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = epa_data_dir
            value = fullfile(misc_emissions_analysis.workspace_dir, 'EPA-Data');
        end
        
        function windows = all_time_periods()
            wc = 2006:2013;
            windows = cell(size(wc));
            for i=1:numel(windows)
                windows{i} = (wc(i)-1):(wc(i)+1);
            end
        end
        
        function labels = all_time_period_labels()
            labels = cellfun(@sprintf_ranges, misc_oh_analysis.all_time_periods, 'uniform', false);
        end
        
        function db = rate_db()
            db_file = fullfile(misc_oh_analysis.epa_data_dir, 'voc_rates.sqlite');
            db = sqlite(db_file);
        end
        
        function [annual_files, years] = list_epa_files()
            annual_files = dirff(fullfile(misc_oh_analysis.epa_data_dir, 'annual*.csv'));
            years = regexp({annual_files.name}, '\d{4}(?=\.csv)', 'match', 'once');
            years = cellfun(@str2double, years);
        end
    
        %%%%%%%%%%%%%%%%%%%%
        % Plotting methods %
        %%%%%%%%%%%%%%%%%%%%
        
        function plot_oh_model_curves(varargin)
            p = advInputParser;
            p.addParameter('locs', misc_emissions_analysis.nine_cities);
            p.addParameter('oh_type', 'invert_hcho');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            oh_type = pout.oh_type;
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.locs);
            
            n_locs = numel(loc_inds);
            n_times = numel(misc_oh_analysis.all_time_periods);
            locs = misc_emissions_analysis.read_locs_file();
            locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            
            % prepare 1 figure per location
            figs = gobjects(n_locs, 1);
            
            all_ax = gobjects(n_locs, 4);
            series = gobjects(n_locs, n_times, 2);
            ax_options = {'xscale', 'log'};
            for i_loc = 1:n_locs
                figs(i_loc) = figure;
                all_ax(i_loc,1) = subplot(2,2,1);
                set(all_ax(i_loc,1), ax_options{:});
                title(sprintf('%s weekdays', locs(i_loc).ShortName));
                xlabel('[NO_x] (ppb)'); ylabel('[OH] (ppt)');
                
                all_ax(i_loc,2) = subplot(2,2,2);
                set(all_ax(i_loc,2), ax_options{:});
                title(sprintf('%s weekends', locs(i_loc).ShortName));
                xlabel('[NO_x] (ppb)');
                
                % for P(HOx), VOCr, and alpha
                all_ax(i_loc,3) = subplot(2,2,3);
                all_ax(i_loc,4) = subplot(2,2,4);
                
                % resize to make top plots bigger
                all_ax(i_loc,1).Position(2) = 0.5;
                all_ax(i_loc,1).Position(4) = 0.4;
                all_ax(i_loc,2).Position(2) = 0.5;
                all_ax(i_loc,2).Position(4) = 0.4;
                all_ax(i_loc,3).Position(4) = 0.25;
                all_ax(i_loc,4).Position(4) = 0.25;
                
                subplot_stretch(1.5,2,'figh',figs(i_loc));
            end
            
            % Load each time period's OH file; get the VOCR, PHOx, alpha,
            % [OH] and [NOx]. Use the first three to make the OH vs NOx
            % steady state curve for that year, use [OH] and [NOx] to mark
            % on the plot where the analysis says that city is.
            year_labels = misc_oh_analysis.all_time_period_labels;
            phox_all = nan(n_locs, n_times, 2);
            vocr_all = nan(n_locs, n_times, 2);
            alpha_all = nan(n_locs, n_times, 2);
            
            for i_yr = 1:n_times
                fprintf('Working on %s\n', year_labels{i_yr});
                years = misc_oh_analysis.all_time_periods{i_yr};
                OH = load(misc_emissions_analysis.oh_file_name(years));
                year_color = misc_oh_analysis.get_time_period_color(years);
                for i_loc = 1:n_locs
                    ii_loc = loc_inds(i_loc);
                    wkday_ax = all_ax(i_loc, 1);
                    wkday_oh = OH.locs_wkday(ii_loc).OH.(oh_type);
                    [series(i_loc, i_yr, 1), phox_all(i_loc, i_yr, 1), alpha_all(i_loc, i_yr, 1), vocr_all(i_loc, i_yr, 1)]...
                        = plot_oh_curve(wkday_ax, wkday_oh, year_color);
                    
                    wkend_ax = all_ax(i_loc, 2);
                    wkend_oh = OH.locs_wkend(ii_loc).OH.(oh_type);
                    [series(i_loc, i_yr, 2), phox_all(i_loc, i_yr, 2), alpha_all(i_loc, i_yr, 2), vocr_all(i_loc, i_yr, 2)]...
                        = plot_oh_curve(wkend_ax, wkend_oh, year_color);
                end
            end
            
            % Go back and add the legend to the top two plots and plot the
            % model parameters (PHOx, alpha, VOCr) on the bottom plots
            for i_loc = 1:n_locs
                legend(all_ax(i_loc,1), series(i_loc, :, 1)', year_labels);
                legend(all_ax(i_loc,2), series(i_loc, :, 2)', year_labels);
                
                plot_model_params(all_ax(i_loc, 3), squeeze(phox_all(i_loc, :, 1)), squeeze(alpha_all(i_loc, :, 1)), squeeze(vocr_all(i_loc, :, 1)));
                plot_model_params(all_ax(i_loc, 4), squeeze(phox_all(i_loc, :, 2)), squeeze(alpha_all(i_loc, :, 2)), squeeze(vocr_all(i_loc, :, 2)));
            end
            
            
            
            function [l, phox, alpha, vocr] = plot_oh_curve(ax, loc_oh, color)
                line_opts = {'color', color, 'linewidth', 2};
                bad_vals = false;
                try
                    nox = loc_oh.nox;
                    oh = loc_oh.oh;
                    
                    phox = loc_oh.phox;
                    alpha = loc_oh.alpha;
                    vocr = loc_oh.vocr;
                catch err
                    if strcmp(err.identifier, 'MATLAB:nonExistentField')
                        fprintf('Have to skip this location/year: %s\n', err.message);
                        bad_vals = true;
                    else
                        rethrow(err)
                    end
                end
                
                if bad_vals || isnan(phox) || isnan(alpha) || isnan(vocr)
                    % create a dummy line to show up in the legend
                    l=line(ax, nan, nan, line_opts{:});
                    phox = nan;
                    alpha = nan;
                    vocr = nan;
                    return
                end
                
                x_nox = logspace(9,12,20);
                y_oh = nan(size(x_nox));
                fprintf('  solving HOx steady state: ');
                parfor i = 1:numel(x_nox)
                    fprintf('*');
                    y_oh(i) = hox_ss_solver(x_nox(i), phox, vocr, alpha);
                end
                fprintf('\n');
                
                l=line(ax, x_nox/2e10, y_oh/2e7, line_opts{:});
                line(ax, nox/2e10, oh/2e7, line_opts{:}, 'marker', 'o', 'markersize', 8, 'linestyle', 'none');
            end
            
            function plot_model_params(ax, loc_phox, loc_alpha, loc_vocr)
                x_years = cellfun(@mean, misc_oh_analysis.all_time_periods);
                x_labels = misc_oh_analysis.all_time_period_labels;
                
                l = gobjects(3,1);
                l(1) = line(ax, x_years, loc_phox / 1e7, 'color', 'k', 'linewidth', 2, 'marker', 'o', 'markersize', 8);
                l(2) = line(ax, x_years, loc_alpha * 100, 'color', 'k', 'linewidth', 2, 'marker', '^', 'markersize', 8);
                l(3) = line(ax, x_years, loc_vocr, 'color', 'k', 'linewidth', 2, 'marker', 'p', 'markersize', 8);
                
                set(ax, 'xtick', x_years, 'xticklabels', x_labels, 'xlim', [min(x_years)-1, max(x_years)+1], 'xticklabelrotation', 45);
                legend(ax, l, {'P(HO_x) (ppt/s)', '\alpha (%)', 'VOC_R (s^{-1})'});
            end
        end
        
        function plot_epa_nox(cities)
            
            [annual_files, years] = misc_oh_analysis.list_epa_files();
            nox = cell(numel(years), numel(cities));
            for i_yr = 1:numel(annual_files)
                tab = readtable(annual_files(i_yr).name);
                for i_city = 1:numel(cities)
                    xx = misc_oh_analysis.find_epa_sites_for_city(cities{i_city}, tab);
                    pp = strcmpi(tab.ParameterName, 'Oxides of nitrogen (NOx)');
                    
                    nox{i_yr, i_city} = misc_oh_analysis.to_double(tab{xx&pp, {'ArithmeticMean'}});
                end
            end
            
            
            for i_city = 1:numel(cities)
                proto_mat1 = nan(numel(years)-2, 1);
                proto_mat2 = nan(numel(years)-2, 2);
                
                nox_means = proto_mat1;
                nox_medians = proto_mat1;
                nox_quants = proto_mat2;
                for i_yr = 1:size(nox_means,1)
                    yy = i_yr:(i_yr+2);
                    this_nox = veccat(nox{yy, i_city});
                    nox_means(i_yr) = nanmean(this_nox);
                    nox_medians(i_yr) = nanmedian(this_nox);
                    nox_quants(i_yr, :) = quantile(this_nox, [0.05, 0.95]);
                end
                
                l = gobjects(3,1);
                figure;
                [~,l(3)] = plot_error_envelope_y(years(2:end-1), nox_quants(:,1), nox_quants(:,2), [1 0.5 0]);
                l(2) = line(years(2:end-1), nox_medians, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
                l(1) = line(years(2:end-1), nox_means, 'color', 'k', 'linestyle', 'none', 'marker', 'p');
                legend(l, {'Mean', 'Median', '5th/95th percentile'});
                ylabel('[NO_x] (ppb)');
                title(cities{i_city});
            end
            
        end
        
        function plot_epa_voc(cities, varargin)
            p = advInputParser;
            p.addParameter('plot_speciated', true);
            p.parse(varargin{:});
            pout = p.Results;
            plot_speciated = pout.plot_speciated;
            
            [annual_files, years] = misc_oh_analysis.list_epa_files();
            
            proto_mat = nan(numel(years), numel(cities));
            proto_cell = cell(numel(years), numel(cities));
            nmocs = proto_cell;
            vocrs = proto_mat;
            vocr_seoms = proto_mat;
            spec_vocrs = proto_cell;
            
            for i_yr = 1:numel(annual_files)
                tab = readtable(annual_files(i_yr).name);
                for i_city = 1:numel(cities)
                    city_name = cities{i_city};
                    xx = misc_oh_analysis.find_epa_sites_for_city(city_name, tab);
                    pp = strcmp(tab.ParameterName, 'Total NMOC (non-methane organic compound)');
                    this_nmoc = tab{xx&pp, {'ArithmeticMean'}};
                    if iscell(this_nmoc)
                        this_nmoc = cellfun(@str2double, this_nmoc);
                    end
                    nmocs{i_yr, i_city} = this_nmoc;
                    [vocrs(i_yr, i_city), vocr_seoms(i_yr, i_city), spec_vocrs{i_yr, i_city}] = misc_oh_analysis.compute_vocr(tab, city_name);
                end
            end
            
            proto_mat1 = nan(numel(years)-2, numel(cities));
            proto_mat2 = nan(numel(years)-2, numel(cities), 2);
            proto_cell = cell(numel(years)-2, numel(cities));
            nmoc_means = proto_mat1;
            nmoc_errs = proto_mat1;
            nmoc_medians = proto_mat1;
            nmoc_quants = proto_mat2;
            
            vocr_means = proto_mat1;
            vocr_errs = proto_mat1;
            vocr_spec = proto_cell;
            vocr_spec_errs = proto_mat1;
            
            avail_spc = cell(1, numel(cities));
            for i_city = 1:numel(cities)
                avail_spc{i_city} = get_avail_species(spec_vocrs(:, i_city));
            end
            
            for i_yr = 1:size(nmoc_means, 1)
                for i_city = 1:numel(cities)
                    yy = i_yr:(i_yr+2);
                    this_nmoc = veccat(nmocs{yy, i_city});
                    nmoc_means(i_yr, i_city) = nanmean(this_nmoc);
                    nmoc_errs(i_yr, i_city) = nanstd(this_nmoc);
                    nmoc_medians(i_yr, i_city) = nanmedian(this_nmoc);
                    nmoc_quants(i_yr, i_city, :) = quantile(this_nmoc, [0.05, 0.95]);
                    
                    vocr_means(i_yr, i_city) = nanmean(vocrs(yy, i_city));
                    vocr_errs(i_yr, i_city) = sqrt(nanmean(vocr_seoms(yy, i_city).^2));
                    [vocr_spec{i_yr, i_city}, vocr_spec_errs(i_yr, i_city)] = merge_speciated_vocrs(spec_vocrs(yy, i_city), avail_spc{i_city});
                end
            end
            
            for i_city = 1:numel(cities)
                figure;
                title(cities{i_city});
                ax = gca;
                %                 plot_error_envelope_y(years(2:end-1), nmoc_quants(:,1), nmoc_quants(:,2));
                %                 line(years(2:end-1), nmoc_medians, 'color', 'k', 'linewidth', 2);
                %                 line(years(2:end-1), nmoc_means, 'color', 'r', 'marker', 'p');
                %                 ylabel('Total NMOC (ppb C)');
                if plot_speciated
                    spc_mat = cat(1, vocr_spec{:,i_city});
                    total_spc_vocr = sum(spc_mat, 2)';
                    yyaxis(ax, 'left');
                    area(years(2:end-1), spc_mat, 'facecolor', 'flat');
                    colormap(jet);
                    common_opts = {'color', 'k', 'linewidth', 2, 'linestyle', '--'};
                    line(years(2:end-1), total_spc_vocr + vocr_spec_errs(:,i_city)', common_opts{:});
                    line(years(2:end-1), total_spc_vocr - vocr_spec_errs(:,i_city)', common_opts{:});
                    ylabel('VOC_R by species (s^{-1})');
                    
                    yyaxis(ax, 'right');
                    common_opts = {'color', 'r', 'linewidth', 2};
                    line(years(2:end-1), nmoc_means(:,i_city), common_opts{:});
                    line(years(2:end-1), nmoc_means(:,i_city) + nmoc_errs(:,i_city), common_opts{:}, 'linestyle', '--');
                    line(years(2:end-1), nmoc_means(:,i_city) - nmoc_errs(:,i_city), common_opts{:}, 'linestyle', '--');
                    ylabel('Total NMOC (ppbC)');
                else
                    yyaxis(ax, 'left')
                    plot_error_envelope_y(years(2:end-1), vocr_means(:,i_city) - vocr_errs(:,i_city),...
                        vocr_means(:,i_city) + vocr_errs(:,i_city), 'r');
                    
                    yyaxis(ax, 'right')
                    plot_error_envelope_y(years(2:end-1), nmoc_quants(:,1), nmoc_quants(:,2), 'r');
                    
                    yyaxis(ax, 'left')
                    line(years(2:end-1), vocr_means, 'color', 'k', 'linewidth', 2);
                    yylabel('VOC_R (s^{-1})');
                    
                    yyaxis(ax, 'right')
                    
                    line(years(2:end-1), nmoc_medians, 'color', 'r', 'linewidth', 2);
                    line(years(2:end-1), nmoc_means, 'color', 'r', 'marker', 'p', 'linestyle', 'none');
                    ylabel('Total NMOC (ppb C)');
                end
            end
            
            function fns = get_avail_species(spec_vocrs)
                % Since its possible that different years may have
                % different measurements, this will create a struct that
                % has only measurements that are in all time periods
                first = 1;
                while isempty(spec_vocrs{first}.Species)
                    first = first + 1;
                    if first > length(spec_vocrs)
                        fns = {};
                        return
                    end
                end
                fns = spec_vocrs{first}.Species;
                for i=(first+1):numel(spec_vocrs)
                    fns = intersect(fns, spec_vocrs{i}.Species);
                end
            end
            
            function [spc, spc_total_err] = merge_speciated_vocrs(spec_vocrs, fns)
                % Since its possible that different years may have
                % different measurements, this will create a struct that
                % has only measurements that are in all time periods
                
                spc = zeros(1, numel(fns));
                spc_total_err = 0;
                for j=1:numel(fns)
                    spc_vocr = 0;
                    for i=1:numel(spec_vocrs)
                        xx_spc = strcmp(spec_vocrs{i}.Species, fns{j});
                        this_vocr = spec_vocrs{i}.VOCR(xx_spc);
                        this_vocr_seom = spec_vocrs{i}.VOCR_SEOM(xx_spc);
                        if ~isnan(this_vocr)
                            spc_vocr = spc_vocr + this_vocr;
                            spc_total_err = spc_total_err + this_vocr_seom.^2;
                        end
                    end
                    spc(j) = spc_vocr;
                end
                spc_total_err = sqrt(spc_total_err);
            end
            
            
        end
        
        function plot_hcho_refsec_diff(year, dow)
            old_file = misc_emissions_analysis.avg_file_name(year, dow, 'hcho');
            [dirname, basename, ext] = fileparts(old_file);
            new_file = fullfile(dirname, 'NewHCHO', [basename, ext]);
            old_vcds = load(old_file);
            old_vcds = old_vcds.daily;
            new_vcds = load(new_file);
            new_vcds = new_vcds.daily;
            
            tmp = veccat(old_vcds.hcho(:), new_vcds.hcho(:));
            cbrange = calc_plot_limits(tmp(~isoutlier(tmp)));
            tmp = veccat(old_vcds.stddev(:), new_vcds.stddev(:));
            cbrange_sd = calc_plot_limits(tmp(~isoutlier(tmp)));
            
            figure;
            subplot_stretch(2,3);
            
            subplot(2,3,1);
            plot_helper(old_vcds.hcho, cbrange, sprintf('Old HCHO Columns (%d)', year));
            
            subplot(2,3,4);
            plot_helper(old_vcds.stddev, cbrange_sd, 'Std. Dev.');
            
            subplot(2,3,2);
            plot_helper(new_vcds.hcho, cbrange, sprintf('New HCHO Columns (%d)', year));

            subplot(2,3,5);
            plot_helper(new_vcds.hcho, cbrange_sd, 'Std. Dev.');
            
            subplot(2,3,3)
            coldel = new_vcds.hcho - old_vcds.hcho;
            plot_helper(coldel, calc_plot_limits(coldel(:), 'diff'), 'Absolute difference', blue_red_cmap);
            
            subplot(2,3,6);
            colrdel = reldiff(new_vcds.hcho, old_vcds.hcho)*100;
            plot_helper(colrdel, calc_plot_limits(colrdel(~isoutlier(colrdel)), 'diff'), 'Percent difference', blue_red_cmap);
            
            function plot_helper(vcds, cblimits, titlestr, cmap)
                if nargin < 4
                    cmap = parula;
                end
                pcolor(old_vcds.lon, old_vcds.lat, vcds); 
                shading flat
                colorbar;
                colormap(gca, cmap);
                caxis(cblimits);
                title(titlestr);
                state_outlines('k');
            end
        end
        
        function fig = plot_covariance_with_oh(varargin)
            p = advInputParser;
            p.addParameter('loc_inds', 1:71);
            p.addParameter('var1', 'oh');
            p.addParameter('var2', 'vocr'); % any field in the OH type structures
            p.addParameter('corr', true); % true = correlation; false = covariance
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_inds = pout.loc_inds;
            main_var = pout.var1;
            co_var = pout.var2;
            do_correlation = pout.corr;
            if do_correlation
                cov_fxn = @corrcoef;
                y_var = 'correlation';
            else
                cov_fxn = @cov;
                y_var = 'covariance';
            end
            
            years = 2006:2013;
            
            oh_types = {'invert', 'invert_hcho', 'invert_hcho_wkday_wkend'};
            type_styles = {'ko', 'b^', 'r*'};
            type_legends = {'Lifetime+SS', 'Lifetime+HCHO SS', 'HCHO SS (VOCR same wkday/wkend)'};
            data = misc_oh_analysis.load_oh_data_for_years(years, 'loc_ids', loc_inds, 'oh_types', oh_types, 'variables', {main_var, co_var});
            var1 = data.(main_var);
            var2 = data.(co_var);
            
            cov_values = nan(numel(loc_inds), numel(oh_types), 2);
            
            for i_loc = 1:numel(loc_inds)
                for i_type = 1:numel(oh_types)
                    for i_dow = 1:2
                        this_oh = var1(:,i_loc,i_type,i_dow);
                        this_var2 = var2(:,i_loc,i_type,i_dow);
                        xx = ~isnan(this_oh) & ~isnan(this_var2);
                        if sum(xx)/numel(this_oh) > 0.5
                            % Only calculate a correlation if we have at
                            % least half of the years
                            tmp_cov = cov_fxn(this_oh(xx), this_var2(xx));
                            cov_values(i_loc, i_type, i_dow) = tmp_cov(1,2); % want the off diagonal term
                        end
                    end 
                end
            end
            
            fig = figure;
            ax1 = subplot(2,1,1);
            hold on
            ax2 = subplot(2,1,2);
            hold on
            
            x = 1:numel(loc_inds);
            for i_type = 1:numel(oh_types)
                plot(ax1, x, cov_values(:, i_type, 1), type_styles{i_type});
                plot(ax2, x, cov_values(:, i_type, 2), type_styles{i_type});
            end
            
            locs = misc_emissions_analysis.read_locs_file();
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(loc_inds);
            locs = misc_emissions_analysis.cutdown_locs_by_index(locs, loc_inds);
            
            legend(ax1, type_legends);
            title(ax1, 'Weekdays');
            legend(ax2, type_legends);
            title(ax2, 'Weekends');
            ax_opts = {'XTickLabel', {locs.ShortName}, 'XTickLabelRotation', 90, 'XLim', [0 x(end)+1], 'YLim', [-1.2, 1.2], 'XTick', x, 'YTick', -1:0.5:1, 'YGrid', 'on'};
            set(ax1, ax_opts{:});
            set(ax2, ax_opts{:});
            
            ylabel_str = sprintf('%s-%s %s', upper(main_var), upper(co_var), y_var);
            ylabel(ax1, ylabel_str);
            ylabel(ax2, ylabel_str);
            subplot_stretch(2,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot helper functions %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function col = get_time_period_color(time_period)
            time_period = mean(time_period);
            all_tps = cellfun(@mean, misc_oh_analysis.all_time_periods);
            col = map2colmap(time_period, min(all_tps), max(all_tps), jet);
        end
        
        function [vocr, vocr_sigma, speciated_vocr] = compute_vocr(epa_table, city)
            % COMPUTE_VOCR Calculate the VOCR from EPA insitu measurements
            %
            % [VOCR, VOCR_SIGMA] = COMPUTE_VOCR(EPA_TABLE, CITY) computes
            % the VOCR and its uncertainty from the table EPA_TABLE
            % (created by loading one of the
            % "annual_conc_by_monitor_yyyy.csv" files with READTABLE) for
            % the city CITY.
            %
            % VOCR is the sum of the average VOCRs for each measured
            % species averaged across all sites in CITY. VOCR_SIGMA is the
            % quadrature sum of the standard errors of the mean for each
            % species across all measurement sites.
            
            % For each species listed in the database, find it in the
            % table, then figure out if its in units of ppbC or ng/m^3. The
            % former is defined by
            % https://www3.epa.gov/ttnamti1/files/ambient/npapsop/sop017.pdf
            % as ppbV * number of carbons (pg. 11). So to get ppb on the
            % way to number density, we need to divide by the number of
            % carbons. We'll handle that by counting carbons in the SMILES
            % string. For the others, we just need the molecular weight to
            % convert. But there's no point in doing any of this if we
            % didn't find an MCM rate constant, so check that first.
            
            db = misc_oh_analysis.rate_db;
            tmp_data = db.fetch('select Parameter, SMILES, MCMRate, MolecularWeight from rates');
            species = tmp_data(:,1);
            smiles = tmp_data(:,2);
            rates = tmp_data(:,3);
            mol_wt = tmp_data(:,4);
            
            vocr = 0;
            vocr_sigma = 0;
            speciated_vocr = cell(0,3);
            missing_species = {};
            
            ll = misc_oh_analysis.find_epa_sites_for_city(city, epa_table);
            
            for i_sp = 1:numel(species)
                if length(rates{i_sp}) <= 1
                    % NULL values for mising rates show up as empty
                    % strings; also some of them I had marked with a "!" to
                    % go back and double check, this skips both.
                    continue
                end
                
                oh_rxn_rate = eval_rate(rates{i_sp});
                
                % Now get the concentrations and units for the current 
                % parameter
                xx = strcmp(epa_table.ParameterName, species{i_sp}) & ll;
                if sum(xx) == 0
                    missing_species{end+1} = species{i_sp};
                    continue
                end
                
                % Double check that we didn't accidentally get sites in two
                % very different cities. Use 2 degree as the cutoff b/c
                % that's how wide the EMG boxes are.
                max_separation = calc_max_dist(misc_oh_analysis.to_double(epa_table.Longitude(xx)),...
                    misc_oh_analysis.to_double(epa_table.Latitude(xx)));
                if max_separation > 2
                    error('Sites chosen for "%s" have a maximum separation of %f degrees - exceeds the maximum allowed', city, max_separation)
                end
                
                species_concs = misc_oh_analysis.to_double(epa_table.ArithmeticMean(xx));
                species_units = unique(epa_table.UnitsOfMeasure(xx));
                if numel(species_units) ~= 1
                    error('Multiple units found (%s) for species "%s" in city "%s"', strjoin(species_units, ', '), species{i_sp}, city);
                else
                    species_units = species_units{1};
                end
                
                switch species_units
                    case 'Parts per billion Carbon'
                        n_carbon = count_carbons(smiles{i_sp});
                        % convert ppbC -> ppb -> molec/cm^3
                        species_concs = species_concs ./ n_carbon .* 1e-9 .* 2e19;
                    case 'Nanograms/cubic meter (25 C)'
                        mw = mol_wt{i_sp};
                        % ng / m^3 -> g / m^3 -> mol/m^3 -> molec./m^3 ->
                        % molec/cm^3
                        species_concs = species_concs .* 1e-9 ./ mw .* 6.022e23 ./ (100^3);
                    otherwise
                        error('No conversion defined for units "%s"', species_units);
                end
                
                % Finally VOCR is the concentration of the VOC in number
                % density times its OH rate constant. Since some species
                % may not be measured at all sites, we'll average each
                % species and compute its standard error of the mean. We'll
                % add the errors in quadrature to get the total uncertainty
                % in the VOCR. Using SEOM to accound for the fact that more
                % sites measuring a particular species reduce the
                % uncertainty of that species.
                
                this_vocr = species_concs .* oh_rxn_rate;
                this_vocr_seom = nanstd(this_vocr) ./ sum(~isnan(this_vocr));
                this_vocr = nanmean(this_vocr);
                
                if ~isnan(this_vocr)  % just in case all concentrations were NaN
                    vocr = vocr + this_vocr;
                    vocr_sigma = vocr_sigma + this_vocr_seom.^2;
                    speciated_vocr = cat(1, speciated_vocr, {species{i_sp}, this_vocr, this_vocr_seom});
                end
            end
            
            if isempty(speciated_vocr)
                % If nothing was found, set VOCR to NaN, not 0.
                vocr = nan;
                vocr_sigma = nan;
            end
            
            vocr_sigma = sqrt(vocr_sigma);
            speciated_vocr = cell2table(speciated_vocr, 'VariableNames', {'Species', 'VOCR', 'VOCR_SEOM'});
            if numel(missing_species) > 0
                missing_species_list = strjoin(missing_species, '\n  * ');
                warning('%d of %d species could not be found for %s:\n  * %s', numel(missing_species), numel(species), city, missing_species_list);
            end
            
            function nc = count_carbons(smiles)
                % Look for a C not followed by a lower case letter (e.g. Cl
                % = chlorine, not a carbon). Include the end of the
                % string.
                nc = numel(regexp(smiles, 'C(?=[^a-z])|C$'));
            end
            
            function rate = eval_rate(rate_expr)
                rate_expr = strrep(rate_expr, '@', '^');
                rate_expr = strrep(rate_expr, 'EXP', 'exp');
                TEMP = 298;
                M = 2e19;
                
                switch rate_expr
                    case 'KMT15'
                        rate = KMT15(TEMP, M);
                    case 'KMT16'
                        rate = KMT16(TEMP, M);
                    case 'KMT17'
                        rate = KMT17(TEMP, M);
                    otherwise
                        rate = eval(rate_expr);
                end
            end
            
            function r_longest = calc_max_dist(lon, lat)
                r_mat = zeros(numel(lon));
                for i=1:numel(lon)
                    for j=(i+1):numel(lon)
                        r = (lon(i) - lon(j)).^2 + (lat(i) - lat(j)).^2;
                        r_mat(i,j) = r;
                    end
                end
                r_longest = max(sqrt(r_mat(:)));
            end
        end
        
        function data = load_oh_data_for_years(years, varargin)
            p = advInputParser;
            p.addParameter('loc_ids', 1:71);
            p.addParameter('oh_types', {}); % leave empty for all
            p.addParameter('variables', {'oh'}); % cell array of fieldnames in each OH struct. If missing, then will be skipped
            
            p.parse(varargin{:});
            pout = p.Results;
            
            loc_inds = misc_emissions_analysis.convert_input_loc_inds(pout.loc_ids);
            oh_types = pout.oh_types;
            variables = pout.variables;
            
            for i_yr = 1:numel(years)
                yr = years(i_yr);
                yr_window = (yr-1):(yr+1);
                fprintf('Loading %s\n', sprintf_ranges(yr_window));
                OH = load(misc_emissions_analysis.oh_file_name(yr_window));
                OH.locs_wkday = misc_emissions_analysis.cutdown_locs_by_index(OH.locs_wkday, loc_inds);
                OH.locs_wkend = misc_emissions_analysis.cutdown_locs_by_index(OH.locs_wkend, loc_inds);
                
                if i_yr == 1
                    if isempty(oh_types)
                        oh_types = fieldnames(OH.locs_wkday(1).OH);
                    end
                    % each array will be years x locs x oh type x days of
                    % week
                    data = make_empty_struct_from_cell(variables, nan(numel(years), numel(loc_inds), numel(oh_types), 2));
                end
                
                for i_loc = 1:numel(loc_inds)
                    for i_type = 1:numel(oh_types)
                        this_type = oh_types{i_type};
                        weekday_oh = OH.locs_wkday(i_loc).OH.(this_type);
                        weekend_oh = OH.locs_wkend(i_loc).OH.(this_type);
                        for i_var = 1:numel(variables)
                            this_var = variables{i_var};
                            if isfield(weekday_oh, this_var)
                                data.(this_var)(i_yr, i_loc, i_type, 1) = weekday_oh.(this_var);
                            end
                            if isfield(weekend_oh, this_var)
                                data.(this_var)(i_yr, i_loc, i_type, 2) = weekend_oh.(this_var);
                            end
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%
        % Misc Methods %
        %%%%%%%%%%%%%%%%
        
        function rows = get_rates()
            fid = fopen(fullfile(misc_oh_analysis.epa_data_dir, 'rate4.csv'));
            tline = fgetl(fid);
            fields = strsplit(tline, '\t');
            tline = fgetl(fid);
            rows = struct([]);
            while ischar(tline)
                tline = regexprep(tline, '\%.+$', '');
                if isempty(tline)
                    tline = fgetl(fid);
                    continue
                end
                line_data = strsplit(tline, '\t');
                line_data = cat(2, line_data, repmat({''}, 1, numel(fields)-numel(line_data)));
                this_line = make_struct_from_field_values(fields, line_data);
                if isempty(rows)
                    rows = this_line;
                else
                    rows(end+1) = this_line;
                end
                tline = fgetl(fid);
            end
            fclose(fid);
        end
        
        function result = get_param_info()
            mcm = misc_oh_analysis.get_rates();
            tab = readtable(fullfile(misc_oh_analysis.epa_data_dir, 'annual_conc_by_monitor_2005.csv'));
            
            params = {'MethodName', 'MetricUsed', 'UnitsOfMeasure'};
            values = struct([]);
            
            for i_var = 1:numel(mcm)
                this_param = mcm(i_var).Parameter;
                xx = strcmp(tab.ParameterName, this_param);
                if sum(xx) == 0
                    fprintf('Could not find %s in the original table\n', this_param);
                    continue
                end
                
                values(end+1).Parameter = this_param;
                
                for p = 1:numel(params)
                    rows = unique(tab{xx, params(p)});
                    if numel(rows) > 1
                        strval = sprintf('(%s)', strjoin(rows, ', '));
                    else
                        strval = rows{1};
                    end
                    values(end).(params{p}) = strval;
                end
            end
            result = struct2table(values);
        end
        
        function xx = find_epa_sites_for_city(city, epa_table)
            loc = misc_oh_analysis.match_loc_to_city(city);
            xx = misc_oh_analysis.find_epa_sites_for_loc(loc, epa_table);
        end
        
        function xx = find_epa_sites_for_loc(loc, epa_table)
            lon = loc.Longitude;
            lat = loc.Latitude;
            radius = mean(loc.BoxSize(3:4));
            
            site_lons = misc_oh_analysis.to_double(epa_table.Longitude);
            site_lats = misc_oh_analysis.to_double(epa_table.Latitude);
            
            xx = (site_lons - lon).^2 + (site_lats - lat).^2 < radius .^2;
        end
        
        function val = to_double(val)
            if iscell(val)
                val = cellfun(@str2double, val);
            end
        end
        
        function loc = match_loc_to_city(city)
            all_locs = misc_emissions_analysis.read_locs_file();
            for i = 1:numel(all_locs)
                if strcmpi(all_locs(i).Location, city) || strcmpi(all_locs(i).ShortName, city)
                    loc = all_locs(i);
                    return
                end
            end
            
            error('Could not find a location matching "%s"', city);
        end
    end
end

