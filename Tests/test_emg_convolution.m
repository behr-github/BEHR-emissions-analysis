function [x, convolved_emg, direct_emg] = test_emg_convolution(f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

if nargin < 1
    A = 2;
    x_0 = 3*3600*5/1000; % x0 = tau * u, use a 3 hr (3*3600 sec) lifetime at 5 m/s winds in km
    mu_x = 20;
    sigma_x = 10;
    B = 1;
    f = [A, x_0, mu_x, sigma_x, B];
elseif ~isvector(f) || ~isnumeric(f) || numel(f) ~= 5
    E.badinput('F must be a 5 element numeric vector, if given');
end

% Define a decently large domain
x = linspace(-100,200,60);

direct_emg = emgfxn_lu(f, x);
conv_fxn = convolved_fit_function(x, gaussian_ld(x, f(4)));
convolved_emg = conv_fxn(f, x);

figure;
plot(x, direct_emg);
hold on
plot(x, convolved_emg);
legend('Direct', 'Convolved');

end

% Copied from fit_line_density on 16 Mar 2018
function e = emgfxn_lu(f,x)
e = f(1)/(2 * f(2)) .* exp( f(3) / f(2) + f(4).^2 / (2*f(2).^2) - x/f(2) )...
    .* erfc( -1/sqrt(2) * ((x-f(3))/f(4) - f(4)/f(2)) ) + f(5);
e(isnan(e)) = Inf;
end

function g = gaussian_ld(x, sigma_x)
g = 1/(sqrt(2*pi)*sigma_x) * exp(-x.^2/(2.*sigma_x.^2));
end