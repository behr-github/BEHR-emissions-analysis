function [ step_lon, step_lat ] = mcmc_choose_center( avg_no2, lon, lat, center_lon, center_lat, radius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
E = JLLErrors;
if ~isequal(size(avg_no2), size(lon)) || ~isequal(size(avg_no2), size(lat))
    E.badinput('AVG_NO2, LON, and LAT must be the same size');
end


% Choose a starting point somewhere within the city radius
[step_lon, step_lat] = choose_start(lon, lat, center_lon, center_lat, radius);

nsteps = 1000;
lon_record = nan(nsteps,1);
lat_record = nan(nsteps,1);
acceptances = false(nsteps,1);
for istep = 1:nsteps
    % Choose the possible next step from a 2D Gaussian distribution
    % centered on the current step
    [next_lon, next_lat] = choose_step(lon, lat, step_lon, step_lat, radius);
    
    % Calculate the acceptance probability, weighting towards higher NO2
    % columns.
    if do_accept(avg_no2, lon, lat, step_lon, step_lat, next_lon, next_lat)
        step_lon = next_lon;
        step_lat = next_lat;
        acceptances(istep) = true;
    end
    lon_record(istep) = step_lon;
    lat_record(istep) = step_lat;
end

% Debugging, control by parameter later
figure; 
pcolor(lon, lat, avg_no2);
shading flat
line(lon_record, lat_record, 'color', 'r');
line(lon_record(1), lat_record(1), 'color', 'r', 'marker', 'x');
line(lon_record(end), lat_record(end), 'color', 'r', 'marker', 'o');
fprintf('%d of %d steps accepted\n', sum(acceptances), nsteps);


end

function [x, y] = choose_start(lon, lat, clon, clat, radius)
r = sqrt((lon(:) - clon).^2 + (lat(:) - clat).^2);
lon = lon(r < radius);
lat = lat(r < radius);
i = randi(numel(lon), 1);
x = lon(i);
y = lat(i);
end

function [x, y] = choose_step(lon, lat, curr_lon, curr_lat, gaussian_sigma)
% Create a 2D gaussian set of weights centered on the current lon/lat, then
% select a point guided by those weights
weights = gauss2d(lon, lat, curr_lon, curr_lat, gaussian_sigma);
isel = randsample(1:numel(lon), 1, true, weights(:));
x = lon(isel);
y = lat(isel);
end

function g = gauss2d(x,y,mu_x,mu_y,s)
% http://mathworld.wolfram.com/GaussianFunction.html
xsqr = (x - mu_x).^2;
ysqr = (y - mu_y).^2;
A = 1./(2*pi*s.^2);
g = A .* exp( -(xsqr + ysqr)./(2*s.^2) );
end

function b = do_accept(no2, lon, lat, x_t, y_t, x_prime, y_prime)
% Prefer to move to higher no2 locations
min_no2 = min(no2(:));
max_no2 = max(no2(:)) - min_no2;
no2_t = interp2(lon, lat, no2, x_t, y_t);
no2_prime = interp2(lon, lat, no2, x_prime, y_prime);

p_t = exp( (no2_t - min_no2)/max_no2 - 1 );
p_prime = exp( (no2_prime - min_no2)/max_no2 - 1 );
% p_t = exp( no2_t - min_no2 );
% p_prime = exp( no2_prime - min_no2 );

b = rand(1) <= p_prime / p_t;
end