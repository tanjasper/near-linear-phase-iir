function [D, dc_inds] = generate_ideal_mag_response(freqs, resp, L, whole)
% [D, freq_inds] = generate_ideal_response(freqs, resp)
%   generates an ideal magnitude only response (zero phase).
%   specify desired responses (from 0 to 1) at freq points freqs (from 0 to
%   1)
%   dc_inds gives the indices of don't care region
%   L specifies number of points for half of unit circle. Function will
%   generally give L+1 points (including pi) or 2L points if whole is set
%   to 'whole'
%   set whole to 'whole' to get response around entire unit circle.
%   otherwise, only half unit circle magnitude will be given

if nargin < 4
    whole = 'half'; % default to generate only for half of unit circle
end
if nargin < 3
    L = 256;
end

D = zeros([L*2 1]);
freq_inds = round(freqs * L) + 1;
dc_inds = 1:freq_inds(1)-1;

for i = 2:length(freqs)
    if resp(i-1) ~= resp(i) % don't care region occurs between two frequencies with different desired responses
        dc_inds = [dc_inds freq_inds(i-1)+1:freq_inds(i)-1];
    end
    D(freq_inds(i-1):freq_inds(i)-1) = resp(i-1);
end

% Finalize by either returning full unit circle or half
if strcmp(whole, 'whole') % full unit circle
    D(end:-1:L+2) = D(2:L); % reflect to other side
    dc_inds = sort([dc_inds -1*dc_inds + 2*L+2]); % reflect to other side
else
    D = D(1:L+1);
end
    
end

