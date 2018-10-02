classdef city_wind_dirs
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static = true)
        function wind_bool = is_wind_good(site, wind_dir)
            % WIND_BOOL = IS_WIND_GOOD( SITE, WIND_DIR ) will return WIND_BOOL,
            % a logical vector the same size as WIND_DIR that indicates whether
            % each wind direction given in WIND_DIR is valid for SITE. SITE 
            % must be a string that specifies a site defined in BAD_WIND_DIRS()
            % attached to this class. WIND_DIR must be a numeric array with
            % values between -180 and +180.

            E = JLLErrors;
            if ~isnumeric(wind_dir) || any(wind_dir(:) < -180 | wind_dir(:) > 180)
                E.badinput('WIND_DIR must be a numeric array with values between -180 and +180');
            elseif ~ischar(site)
                E.badinput('SITE must be a string')
            end

            wind_bool = true(size(wind_dir));
            bad_wind = city_wind_dirs.bad_wind_dirs(site);
            for a=1:numel(wind_dir)
                if any(wind_dir(a) > bad_wind(:,1) & wind_dir(a) < bad_wind(:,2))
                    wind_bool(a) = false;
                end
            end
        end

        function bad_dirs = bad_wind_dirs(site)
            % BAD_WIND_DIRS( SITE ) Given a site name, SITE, as a string,
            % will return a n-by-2 array of wind directions that should not
            % be used for that site, usually because there is a secondary
            % source downwind. Each row in the array will be [min, max] so
            % that wind directions >= min and <= max should be rejected. If
            % there are no bad directions, it will be an empty matrix. These
            % directions are specified in degrees CCW from east; and should
            % be on the range [-180, +180]
            
            E = JLLErrors;
            
            % Switch based on the city name. This should match the short
            % name in trend_locations.nc
            switch site
                case 'Atlanta'
                    bad_dirs = [-112.5, 0];
                otherwise
                    E.badinput('No wind directions defined for SITE == "%s"', site);
            end
            
            E.addCustomError('bad_wind_def', 'The bad_dirs array for %s is defined incorrectly; %s.');
            if size(bad_dirs,2) ~= 0 && size(bad_dirs,2) ~= 2
                E.callCustomError('bad_wind_def', site, 'it is not 2 long in the second dimension');
            elseif any(bad_dirs(:,1) > bad_dirs(:,2))
                E.callCustomError('bad_wind_def', site, 'one or more of the entries in the first column are greater than the corresponding entry in the second');
            elseif any(bad_dirs(:) < -180 | bad_dirs(:) > 180)
                E.callCustomError('bad_wind_def', site, 'one or more of the entries are outside the range [-180, +180]');
            end
        end
    end
    
end

