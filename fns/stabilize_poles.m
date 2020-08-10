function [a] = stabilize_poles(a)
% [a] = stabilize_poles(a)
%   Take denominator polynomial coefficients, take poles, put all poles in
%   unit circle, and return new polynomial coefficients

poles = roots(a);
for k = 1:length(poles)
    if abs(poles(k)) > 1
        poles(k) = 1/poles(k); % if unstable pole, take inverse
    end
end
a = poly(poles);
        
end

