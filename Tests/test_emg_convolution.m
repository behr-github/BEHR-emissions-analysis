function [x, convolved_emg, direct_emg] = test_emg_convolution(f)
%TEST_EMG_CONVOLUTION Test the numerical convolution lifetime fitting approach
%   [ X, CONVOLVED_EMG, DIRECT_EMG ] = TEST_EMG_CONVOLUTION() This function
%   will test the numerical convolution approach from Fei et al. 2016 (doi:
%   10.5194/acp-16-5283-2016) by reproducing an exponentially modified
%   Gaussian with it. This will generate a default set of fitting
%   parameters and then generate an EMG function both with the analytical
%   version from Lu et al. 2015 (doi: 10.5194/acp-15-10367-2015) and by
%   convolving a Gaussian source distribution numerically. It will plot
%   both and try to fit both with fit_line_density() (using the Lu fitting
%   function) and print the original fit parameters used to make the EMG
%   functions, and the fitted parameters to both functions. If all goes
%   well, both should be nearly identical. It will return the
%   x-coordinates, numerically convolved EMG function, and analytical EMG
%   function as X, CONVOLVED_EMG, and DIRECT_EMG, respectively.
%
%   [ X, CONVOLVED_EMG, DIRECT_EMG ] = TEST_EMG_CONVOLUTION( F ) allows you
%   to specify the fit parameters as the five element vector [a, x_0, mu_x,
%   sigma_x, B].

E = JLLErrors;

if nargin < 1
    init_fit.a = 1;%2;
    init_fit.x_0 = 3*3600*5/1000; % x0 = tau * u, use a 3 hr (3*3600 sec) lifetime at 5 m/s winds in km
    init_fit.mu_x = 20;
    init_fit.sigma_x = 10;
    init_fit.B = 0;%1;
    f = struct2array(init_fit);
elseif ~isvector(f) || ~isnumeric(f) || numel(f) ~= 5
    E.badinput('F must be a 5 element numeric vector, if given');
else
    init_fit = make_struct_from_field_values({'a','x_0','mu_x','sigma_x','B'}, f);
end

% Define a decently large domain
x = linspace(-100,200,60);

direct_emg = emgfxn_lu(f, x);
conv_fxn = convolved_fit_function(x, gaussian_ld(x, f(4)));
convolved_emg = conv_fxn(f, x);

% Beirle doesn't have the x_0 in the prefactor - what happens if I take it
% out of the Lu function?

% Test what the fit says the parameters are
[ffit_direct, fit_emg_direct] = fit_line_density(x, direct_emg,'none');
[ffit_conv, fit_emg_conv] = fit_line_density(x, convolved_emg,'none');
%[ffit_conv2, fit_emg_conv2] = fit_line_density(x, convolved_emg, 'none', 'emgtype', 'conv-gaussian');

fprintf('Original fit params: %s\n', struct2string(init_fit));
fprintf('Direct function fit params: %s\n', struct2string(ffit_direct));
fprintf('Convolved function fit params: %s\n', struct2string(ffit_conv));
%fprintf('Convolved function fit using same function in fitting: %s\n', struct2string(ffit_conv2));

figure;
plot(x, direct_emg, 'bo');
hold on
plot(x, convolved_emg, 'ro');
plot(x, fit_emg_direct, 'b--');
plot(x, fit_emg_conv, 'r--');
%plot(x, fit_emg_conv2, 'r-.');
legend('Orig. direct', 'Orig. convolved','Fit direct','Fit convolved')%,'Fit convolved with convolved');

end

% Copied from fit_line_density on 16 Mar 2018
function e = emgfxn_lu(f,x)
e = f(1)/2 .* exp( f(3) / f(2) + f(4).^2 / (2*f(2).^2) - x/f(2) )...
    .* erfc( -1/sqrt(2) * ((x-f(3))/f(4) - f(4)/f(2)) ) + f(5);
e(isnan(e)) = Inf;
end

function g = gaussian_ld(x, sigma_x)
g = 1/(sqrt(2*pi)*sigma_x) * exp(-x.^2/(2.*sigma_x.^2));
end