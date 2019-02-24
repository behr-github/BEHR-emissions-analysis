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
p.addParameter('variables', {'no2_vcds'});
p.addParameter('avg_levels', 'all');
p.parse(varargin{:});
pout = p.Results;

variables = pout.variables;
avg_levels = pout.avg_levels;

avg_proc.no2_vcds = struct('variables', {{'no2_ndens', 'z', 'z_center'}}, 'proc_fxn', @wrf_no2_vcd);

dvec = make_datevec(start_dates, end_dates);
WRFFiles = BEHRMatchedWRFFiles('region', 'us');
for d=1:numel(dvec)
    fprintf('Working on %s\n', datestr(dvec(d)));
    day_avg = wrf_time_average(dvec(d), dvec(d), variables, 'processing', avg_proc, 'matched_wrf_files', WRFFiles);
    [xloncorn, xlatcorn] = wrf_grid_corners(day_avg.XLONG, day_avg.XLAT);
    Data = struct('Latitude', day_avg.XLONG, 'Longitude', day_avg.XLAT,...
        'FoV75CornerLongitude', xloncorn, 'FoV75CornerLatitude', xlatcorn,...
        'Areaweight', ones(size(day_avg.XLONG)));
    for i_var = 1:numel(variables)
        this_var = variables{i_var};
        value = day_avg.(this_var);
        if ndims(value) > 3
            E.notimplemented('Value averaged from WRF files has >3 dimensions');
        elseif ~ismatrix(value) % if 3d
            if ischar(avg_levels) && strcmpi(avg_levels, 'all')
                levels = 1:size(value,3);
            elseif ischar(avg_levels)
                E.badinput('The only string recognized for "avg_levels" is "all"')
            else
                levels = avg_levels;
            end
            if ~isempty(levels)
                value = nanmean(value(:,:,avg_levels), 3);
            end
        end
        Data.(this_var) = value;
    end
    
    save_name = sprintf('WRF_PseudoBEHR_%04d%02d%02d.mat',year(dvec(d)),month(dvec(d)),day(dvec(d)));
    fprintf('    Saving %s\n',fullfile(save_path,save_name));
    save(fullfile(save_path,save_name),'Data');
end


end

function vcd = wrf_no2_vcd(Wrf)
perm_vec = [3 1 2 4];

z = permute(Wrf.z, perm_vec);
z_center = permute(Wrf.z_center, perm_vec);
no2_ndens = permute(Wrf.no2_ndens, perm_vec);

sz = size(no2_ndens);
% wrf_day_avg needs time in the 4th dimension
vcd = nan([sz(2:3),1,sz(4)]);

for i=1:numel(vcd)
    no2_interp = nan(size(z(:,i)));
    no2_interp(2:end-1) = interp1(z_center(:,i), no2_ndens(:,i), z(2:end-1,i));
    % use constant extrapolation
    no2_interp(1) = no2_ndens(1);
    no2_interp(end) = no2_ndens(end);
    % z in meters, need to integrate in cm
    vcd(i) = trapz(z(:,i)*100, no2_interp);
end


end
