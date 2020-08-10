function [err_val, mean_err_val] = mag_err(b, a, D, dc_inds)
% [err_val] = mag_err(b, a, D, dc_inds)
%   D should only be from 0 to pi

if nargin < 4
    dc_inds = [];
end

% ensure D is a row vector
if size(D,1) == 1
    D = permute(D, [2 1]);
end

L = length(D)-1;

dc_inds(dc_inds > length(D)) = []; % remove dc_inds that go beyond signal

% Obtain freq response from 0 to pi
H = freqz(b,a, 2*L, 'whole');
H = H(1:L+1);

% remove don't care region
H(dc_inds) = [];
D(dc_inds) = [];

err_val = norm(abs(H) - abs(D));
mean_err_val = norm(abs(H) - abs(D)) / length(H);

end

