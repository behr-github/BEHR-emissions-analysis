function c = conv_trapz(x, f, g, method)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('method','var')
    method = 'full';
end

x = x(:);
f = f(:);
g = g(:);

% The definition of a convolution is:
%
%   (f * g)(t) = \int_{-\infty}^{+\infty} f(\tau) g(t - \tau) d\tau
%
% in other words, the valus of the convolution f * g at t is equal to the
% integral of f times g shifted by t and reversed.

% First pad g with enough zeros that we can shift it such that their last
% nonzero values line up with x1 and y1's first values and vice versa. We
% actually flip g to account for it being mathematically reversed by the
% -\tau in the function.
orig_length = length(g);
zero_padding = zeros(orig_length,1);
g = flipud([zero_padding; g; zero_padding]);
g_orig_start = orig_length + 1;

% Now here, c is our convolved function. So as (f*g)(t) is the integral
% when g is shifted by t, c(i) is when the indices of g are shifted by i.
c = nan(size(g));
n = numel(f);
for i=-orig_length:orig_length
    gstart = g_orig_start + i;
    gend = gstart + orig_length - 1;
    cidx = i + orig_length + 1; % transform i so that the index into c starts at 1
    c(cidx) = trapz(x, f .* g(gstart:gend));
end

if strcmpi(method,'same')
    c = c(g_orig_start:(g_orig_start+orig_length-1));
elseif ~strcmpi(method,'full')
    % Full convolution returns c unaltered, as in CONV().
    E.notimplemented('Convolution method "%s" not implemented', method)
end

end

