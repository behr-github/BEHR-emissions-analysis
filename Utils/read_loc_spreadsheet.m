function [ locs ] = read_loc_spreadsheet(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mydir = fileparts(mfilename('fullpath'));
spreadsheet_file = fullfile(mydir, '..', 'Workspaces', 'SiteData', 'trend_locations.xlsx');
[~,sheets] = xlsfinfo(spreadsheet_file);

locs = [];

for a=1:numel(sheets)
    [~,~,raw] = xlsread(spreadsheet_file, sheets{a});
    % Only take columns with headers, some rows have unneeded notes at the
    % end
    xx = cellfun(@(x) ischar(x), raw(1,:)); 
    
    % Stop at a blank line - some sheets have notes at the end
    nans = cellfun(@(x) ~ischar(x) && isnan(x), raw);
    last_line = find(all(nans,2),1) - 1;
    if isempty(last_line)
        last_line = size(raw,1);
    end
    
    header = raw(1,xx);
    raw = raw(2:last_line,xx);
    
    % Convert the box size from a string to a numeric array
    box_size_ind = strcmpi(header, 'BoxSize');
    if sum(box_size_ind) == 0
        error('trend_locations:missing_category', 'The category "BoxSize" is not present in the Excel file')
    else
        header{end+1} = 'SiteType';
        raw(:,box_size_ind) = convert_box_size(raw(:,box_size_ind), sheets{a});
    end
    
    raw(:,end+1) = repmat({sheets{a}},size(raw,1),1);
    these_locs = cell2struct(raw,header,2);
    locs = cat(1, locs, these_locs);
end

end

function cell_out = convert_box_size(cell_in, location_type)
cell_out = cell(size(cell_in));
for a=1:numel(cell_in)
    if isnan(cell_in{a})
        if strcmpi(location_type, 'cities')
            cell_out{a} = [1 2 1 1];
        elseif strcmpi(location_type, 'powerplants')
            cell_out{a} = [0.5 1 0.5 0.5];
        elseif strcmpi(location_type, 'ruralareas')
            cell_out{a} = nan(1,4);
        else
            error('trend_locations:unknown_site_type','No default box size defined for site type "%s"', location_type);
        end
    else
        value = regexprep(cell_in{a},'[^\d\s\.]', '');
        cell_out{a} = str2double(strsplit(value));
    end
end
end