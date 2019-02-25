function [  ] = preproc_WRF2Data( start_dates, end_dates, save_path, varargin )
%PREPROC_WRF2DATA Process WRF BEHR output into Data structures
%   PREPROC_WRF2DATA( START_DATES, END_DATES, SAVE_PATH ) Will process WRF
%   data for every day in the time range(s) defined by START_DATES and
%   END_DATES, which may be any format understood by MAKE_DATEVEC. NO2 VCDs
%   computed from the WRF files that match up with OMI overpass times will
%   be computed and a day average saved as "WRF_PseudoBEHR_yyyymmdd.mat" in
%   the directory specified by SAVE_PATH.
%
%   Parameters:
%
%       'variables' - a cell array of WRF variables to store in the output
%       files. Any variable present in the WRF file or which can be
%       computed by READ_WRF_PREPROC is allowed. In addition 'no2_vcd',
%       which is the only default variable included, will compute NO2 VCDs
%       using INTEGPR2.
%
%       'avg_levels' - controls what levels of 3D WRF variables are
%       averaged to produce 2D output arrays. 'all' (default) will average
%       over all levels, otherwise this should be a vector specifying the
%       levels by index, e.g. 1:5 to average over the first five levels.
%       You may pass an empty array to skip averaging.

E = JLLErrors;

p = advInputParser;
p.addParameter('variables', {'no2_vcd'});
p.addParameter('avg_levels', 'all');

variables = pout.variables;
avg_levels = pout.avg_levels;

data_vars = veccat({'Date', 'Longitude', 'Latitude', 'FoV75CornerLongitude', 'FoV75CornerLatitude', 'Areaweight'}, variables);

avg_proc.no2_vcd = struct('variables', {'no2', 'pres'}, 'proc_fxn', @wrf_no2_vcd);

dvec = make_datevec(start_dates, end_dates);
WRFFiles = BEHRMatchedWRFFiles('region', 'us');
for d=1:numel(dvec)
    fprintf('Working on %s\n', datestr(dvec(d)));
    [todays_wrf_files, behr_data] = WRFFiles.get_files_for_date(dvec(d));
    Data = repmat(make_empty_struct_from_cell(data_vars), size(behr_data));
    for i_orbit = 1:numel(behr_data)
        swath_lon_edge = edge_to_vector(behr_data(i_orbit).Longitude);
        swath_lat_edge = edge_to_vector(behr_data(i_orbit).Latitude);
        wrf_lon = ncread(todays_wrf_files{i_orbit}, 'XLONG');
        wrf_lat = ncread(todays_wrf_files{i_orbit}, 'XLAT');
        [wrf_loncorn, wrf_latcorn] = wrf_grid_corners(wrf_lon, wrf_lat);
        
        xx = inpolygon(wrf_lon, wrf_lat, swath_lon_edge, swath_lat_edge);
        Data(i_orbit).Date = datestr(dvec(d));
        Data(i_orbit).Longitude = wrf_lon(xx);
        Data(i_orbit).Latitude = wrf_lat(xx);
        Data(i_orbit).FoV75CornerLongitude = wrf_loncorn(:,xx);
        Data(i_orbit).FoV75CornerLatitude = wrf_latcorn(:,xx);
        Data(i_orbit).Areaweight = ones(size(Data(i_orbit).Longitude));
        for i_var = 1:numel(variables)
            this_var = variables{i_var};
            if isfield(avg_proc, this_var)
                WRF = read_wrf_vars('', todays_wrf_files{i_orbit}, avg_proc.(this_var).variables, 'squeeze', 'as_struct');
                value = avg_proc.(this_var).proc_fxn(WRF);
            else
                value = read_wrf_preproc(todays_wrf_files{i_orbit}, this_var);
            end
            
            if ndims(value) > 3
                E.notimplemented('WRF array with > 3 dimensions')
            elseif ischar(avg_levels) && strcmpi(avg_levels, 'all')
                value = nanmean(value, 3);
            elseif ischar(avg_levels)
                E.badinput('The only recognized string for "avg_levels" is "all"')
            else
                value = nanmean(value(:,:,avg_levels),3);
            end
            
            Data(i_orbit).(this_var) = value;
        end
    end
    save_name = sprintf('WRF_PseudoBEHR_%04d%02d%02d.mat',year(dvec(d)),month(dvec(d)),day(dvec(d)));
    fprintf('    Saving %s\n',fullfile(save_path,save_name));
    save(fullfile(save_path,save_name),'Data');
end


end

function vcd = wrf_no2_vcd(Wrf)
no2 = permute(Wrf.no2, [3 1 2]);
pres = permute(Wrf.pres, [3 1 2]);
sz = size(no2);
vcd = nan(sz(2:3));
for i=1:numel(no2)
    vcd(i) = integPr2(no2(:,i)*1e-6, pres(:,i), pres(1,i));
end
end