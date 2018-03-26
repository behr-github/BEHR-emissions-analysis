function [ OMI ] = rotate_plume( Data, center_lon, center_lat, theta, varargin )
%ROTATE_PLUME Code to rotate NO2 plumes to align them by wind direction
%   Required inputs:
%       Data - should be the Data structure created by read_omno2 or
%       BEHR_main.m Only pass one top-level element at a time, i.e. if Data
%       is a 1x4 structure, pass in Data(i) where i is 1, 2, 3, or 4.
%
%       center_lon/center_lat - the center coordinates that the rotation
%       should be carried out around. Must be scalar numbers between -180
%       and 180 and -90 and 90 respectively.
%
%       theta - the angle the plume makes relative to a vector due east
%       from the source. That is, a plume being advected to the northwest
%       will have a theta of ~135 degrees.  This should be degrees (not
%       radians). This might be derived based on the wind direction over
%       the source, the major axis of an ellipse fit to the extent of the
%       plume, or another method.
%
%   Optional inputs:
%       rel_box_corners - a 1x4 vector setting the size of the box to
%       contain the plume. The values are: degrees west of center, degrees
%       east of center, degrees south of center, degrees north of center.
%       Defaults to [2 4 2 2], i.e. the box will be 6 deg E/W by 4 deg N/S
%       offset so that it extends 2 deg further E than west.
%
%       vza_crit - maximum viewing zenith angle allowed (in degrees).
%       Default is 60.
%
%       loncorn, latcorn - the field names in Data to use for longitude and
%       latitude corners of the pixels. Defaults are 'FoV75CornerLongitude'
%       and 'FoV75CornerLatitude' respectively.
%
%       pixels_in_box - if true, only returns the logical array (the same
%       size as Data.Longitude and Data.Latitude) that indicates which
%       pixels are inside the box and pass the VZA criterion.


%   Original code from Luke Valin:
%
%         xBOX = [-2 4 4 -2];
%         yBOX = [-2 -2 2  2];
%         V = hypot(windX(date_i),windY(date_i));
%         T = theta(date_i);
%         for corner = 1 : 4
%          out =   [cos(T) -sin(T); sin(T) cos(T) ]*[xBOX(corner); yBOX(corner)];
%          xBp(corner)= out(1)+pp_lon;
%          yBp(corner)= out(2)+pp_lat;
%         end
%
%         v_ind = find(inpolygon(swath_lon,swath_lat,xBp,yBp))
%
%                          if length(v_ind)<=1
%                             vcd_tot(date_i,:,:)=nan(161,241);
%                             continue
%                          end
%
%
%         swath_X = nan(1,length(v_ind));
%         swath_Y = nan(1,length(v_ind));
%         for t_i =1:length(v_ind)
%
%             out = [cos(-T) -sin(-T); sin(-T) cos(-T) ]*[swath_lon( v_ind(t_i))-pp_lon; swath_lat(v_ind(t_i))-pp_lat];
%             swath_X(t_i)=out(1);
%             swath_Y(t_i)=out(2);
%
%         end
%
%   The idea is to define a 4 deg lat x 6 deg lon box that is assumed to
%   capture a plume.  It is defined along the x-axis, with 2 deg to the
%   west of the center point and 4 deg east. This is then rotated to align
%   with the wind direction defined by theta (which should be degrees CCW
%   from east, i.e. a normal definition of theta in polar coordinates). It
%   then identifies which pixels fall in that box and rotates them back to
%   the x-axis.
%
%   There are two potential ways to implement this with BEHR data. The
%   simplest approach would be to rotate the gridded version in the OMI
%   structure and simply assume that rotating the borders of the grid cells
%   is unecessary.  The potential problem with that is that it will not
%   ensure that the resulting grid lies well on the x-y plane, which would
%   make integration across the plume to get a line density more difficult.
%
%   The better solution I think is to rotate the native pixels and then
%   grid them, thus ensuring that the grid is defined along the x-y axes.
%   This will make the following steps easier, but will require some more
%   careful identification of which pixels fall into the box and rotation
%   of the pixel corners as well.
%
%
%   Josh Laughner <joshlaugh5@gmail.com> 4 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addOptional('rel_box_corners',[]);
p.addParameter('vza_crit',60);
p.addParameter('loncorn', 'FoV75CornerLongitude');
p.addParameter('latcorn', 'FoV75CornerLatitude');
p.addParameter('pixels_in_box', false);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;
rel_box_corners = pout.rel_box_corners;
vza_crit = pout.vza_crit;
loncorn_field = pout.loncorn;
latcorn_field = pout.latcorn;
only_pixels_in_box = pout.pixels_in_box;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

E = JLLErrors;

if ~isstruct(Data) || ~isscalar(Data) || any(~ismember({'Longitude','Latitude',loncorn_field,latcorn_field}, fieldnames(Data)))
    E.badinput('Data must be a scalar structure with the fields Longitude, Latitude, %s, and %s', loncorn_field, latcorn_field)
end
if ~isscalar(center_lon) || ~isnumeric(center_lon) || center_lon > 180 || center_lon < -180
    E.badinput('center_lon must be a scalar numeric value between -180 and +180')
elseif ~isscalar(center_lat) || ~isnumeric(center_lat) || center_lat > 180 || center_lat < -180
    E.badinput('center_lat must be a scalar numeric value between -90 and +90')
elseif ~isscalar(theta) || ~isnumeric(theta) || theta > 180 || theta < -180
    E.badinput('theta must be a scalar numeric value between -180 and +180')
end

if any(isnan(rel_box_corners))
    E.badinput('REL_BOX_CORNERS cannot contain NaNs')
elseif ~isempty(rel_box_corners)
    if numel(rel_box_corners) ~= 4 || ~isnumeric(rel_box_corners) || ~isvector(rel_box_corners)
        E.badinput('rel_box_corners (if given) must be a 4 element numeric vector')
    end
    x_box = [-rel_box_corners(1), rel_box_corners(2), rel_box_corners(2), -rel_box_corners(1)];
    y_box = [-rel_box_corners(3), -rel_box_corners(3), rel_box_corners(4), rel_box_corners(4)];
else
    x_box = [-2 4 4 -2];
    y_box = [-2 -2 2 2];
end

if ~isscalar(vza_crit) || ~isnumeric(vza_crit)
    E.badinput('vza_crit must be a numeric scalar')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% First rotate the box to align along theta and find any pixels that
% overlap it at all.
R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
x_box_rot = nan(size(x_box));
y_box_rot = nan(size(y_box));
for corner=1:4
    out = R * [x_box(corner); y_box(corner)];
    x_box_rot(corner) = out(1) + center_lon;
    y_box_rot(corner) = out(2) + center_lat;
end

% Check if this orbit has any pixels within the target box, if not, return
% a dummy structure indicating so
in_check = inpolygon(Data.Longitude, Data.Latitude, x_box_rot, y_box_rot);

% If there are no pixels in the box, then we can exit now. Format the
% output appropriately whether the logical array or the gridded data
% structure was requested.
if ~any(in_check(:))
    if only_pixels_in_box
        OMI = in_check;
    else
        
        OMI.Longitude = [];
        OMI.Latitude = [];
    end
    return;
end


% Remove pixels with VZA greater than the specified criteria which defaults
% to 60 degrees. Do any data field that is the 2D shape, and avoid any flag
% fields, because they will cause a Python error on the CVM side.
vv = Data.ViewingZenithAngle > vza_crit;
fns = fieldnames(Data);
for f=1:numel(fns)
    if isequal(size(Data.(fns{f})), size(vv)) && isnumeric(Data.(fns{f})) && ~ismember(fns{f}, BEHR_publishing_gridded_fields.flag_vars)
        Data.(fns{f})(vv) = nan;
    end
end

% If requested to just return which pixels will be used, combine the test
% for which pixels are in the box and which pass the VZA crit test, then
% return before doing the gridding.
if only_pixels_in_box
    OMI = in_check & ~vv;
    return
end


% Rotate the pixels back to the x-axis (only need to change the lon/lat fields)
R = [cosd(-theta), -sind(-theta); sind(-theta), cosd(-theta)];
for a=1:numel(Data.Longitude)
    out = R * [Data.Longitude(a) - center_lon; Data.Latitude(a) - center_lat];
    Data.Longitude(a) = out(1) + center_lon;
    Data.Latitude(a) = out(2) + center_lat;
    for b=1:4
        out = R * [Data.(loncorn_field)(b,a) - center_lon; Data.(latcorn_field)(b,a) - center_lat];
        Data.(loncorn_field)(b,a) = out(1) + center_lon;
        Data.(latcorn_field)(b,a) = out(2) + center_lat;
    end
end


% Finally grid the data to a 0.05 x 0.05 degree grid.
lonmin = center_lon + x_box(1);  lonmax = center_lon + x_box(3);
latmin = center_lat + y_box(1);  latmax = center_lat + y_box(3);
resolution = 0.05;
BoxGrid = GlobeGrid(resolution, 'domain', [lonmin, lonmax, latmin, latmax]);

OMI = psm_wrapper(Data, BoxGrid, 'only_cvm', true, 'DEBUG_LEVEL', DEBUG_LEVEL);


end

