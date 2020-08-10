function [rs, grp_delay] = lin_phase_r2(b, a, L, pb_inds, showfig)
% [err_val] = lin_phase_err(b, a, pb_inds)
%   Finds the l2 error from the best fitting line to the phase given by b,a
%   pb_inds is the indices of the passband (leave empty to perform for
%   entire response)
%   Generally, we don't care about phase in stop band (because magnitude is
%   0)

if nargin < 5
    showfig = false;
end
if nargin < 3
    L = 512;
end
if nargin < 4
    pb_inds = 1:L;
end

H = freqz(b,a,L);
ph = unwrap(angle(H));
line_x = 0:L-1;

% check passband only
ph = ph(pb_inds);
line_x = line_x(pb_inds)';

% slope of best fitting line
m = line_x \ ph;
grp_delay = -1 * (m/pi) * L;
% err_val = norm(ph - line_x * m);

% find r-squared value
% ref: https://onlinecourses.science.psu.edu/stat501/node/255/
ssr = sum((ph - mean(ph)).^2);
sse = sum((ph - line_x*m).^2);
rs = ssr/(ssr+sse);

if showfig
    figure
    plot(unwrap(angle(H)));
    hold on, plot(line_x * m);
    title(sprintf('R^2: %f%%', rs * 100));
    legend('System phase', 'best fit line');
end

end

