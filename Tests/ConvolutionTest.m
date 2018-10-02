classdef ConvolutionTest < matlab.unittest.TestCase
    % CONVOLUTIONTEST - Tests my convolution for EMG emissions analysis
    %   Run with run(ConvolutionTest())
    
    methods(Test)
        function testBoxcarConvolution(testCase)
            % The convolution of two boxcar functions should have a maximum
            % when the two fully overlap, which if we test by convolving
            % the same boxcar with itself, should be 0. Also the max itself
            % should be the result of integrating the function times
            % itself.
            x = -2:0.01:2;
            y = 2*double(abs(x) <= 1);
            
            cprime = conv_trapz(x,y,y);
            [max_val, max_idx] = max(cprime);
            testCase.verifyEqual(max_val, trapz(x,y.*y));
            testCase.verifyEqual(x(max_idx), 0);
        end
        
        function testGaussianConvolution(testCase)
            % From http://mathworld.wolfram.com/Convolution.html, the
            % convolution of two Gaussians is a Gaussian with a variance
            % equal to the sum of the two convolved Gaussian's variances.
            % So test this by convolving the a Gaussian with itself and
            % compare it to one that has twice the variance (sigma.^2).
            x = -2:0.01:2;
            % Define a gaussian by its x coordinates and variance
            % (sigma.^2).
            g = @(x,s2)  1./(sqrt(2*pi*s2))*exp(-x.^2 / (2*s2));
            sig = 0.0625;
            cprime = conv_trapz(x, g(x,sig), g(x,sig));
            
            % The two vectors will not be exactly equal due to floating
            % point error. I picked 10^-10 as the criteria just as it is
            % arbitrarily small.
            residual = sqrt(sum((cprime(:) - g(x(:),2*sig)).^2));
            testCase.verifyLessThanOrEqual(residual, 1e-10);
        end
    end
end

