function cov_mat = emg_cov_mat(hessian_mat, no2_x, no2_ld, emgfit)
%EMG_COV_MAT Compute the covariance matrix for the EMG fitting
%   COV_MAT = EMG_COV_MAT( HESSIAN_MAT, NO2_X, NO2_LD, EMGFIT ) Computes
%   the covariance matrix for an EMG fit. Requires the Hessian matrix
%   returned by the fitting process, the x-coordinates for the fit and the
%   line density (NO2_X), the line densities (NO2_LD), and the vector
%   describing the fit (EMGFIT).
%
%   The covariance matrix is computed as the inverse of the Hessian matrix
%   returned multiplied by the sum of squared errors. The inverse of the
%   Hessian is an approximation to the covariance matrix (Dovi 1991, Appl.
%   Math. Lett. Vol 4, No 1, pp. 87-90); for a single-equation system, it
%   should also be multiplied by the standard squared error of the fit
%   (basically the sum of squared residuals divided by the number of
%   degrees of freedom, see also Ch. 7 of "Nonlinear Parameter Estimation"
%   by Yonathan Bard).

n_ld = sum(~isnan(no2_x) & ~isnan(no2_ld));
sse = nansum2((emgfit - no2_ld).^2) * 1/(n_ld - 5);
cov_mat = sse * inv(hessian_mat);
end

