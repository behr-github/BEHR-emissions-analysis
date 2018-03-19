function emgfxn = convolved_fit_function(x_slow, slow_ld)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

emgfxn = @convolved_fxn;

    function emg = convolved_fxn(f,x)
        if ~isequal(x, x_slow)
            error('convolved_emg_fit:x_mismatch', 'New x coordinates do not match those passed for the slow line densities')
        end
        
        a = f(1);
        x_0 = f(2);
        mu_x = f(3);
        sigma_x = f(4);
        B = f(5);
        
        exp_component = zeros(size(x));
        exp_component(x > mu_x) = exp( -(x(x>mu_x)-mu_x)/x_0 );
        emg = conv(exp_component, a .* slow_ld, 'same') + B;
    end

end

