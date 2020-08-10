function [b, a, best_err] = gauss_newton_iir(D, b_init, a_init, dc_inds, num_iters)
% Generates an IIR filter based on the Gauss-Newton method
%   [b,a,best_err] = gauss_newton_iir(D,b_init,a_init,dc_inds,num_iters)
%   returns a filter whose response approximates the complex response D as
%   closely as possible. The method is initialize with b_init and a_init.
%   dc_inds specifies don't-care indices of D, and num_iters specifies the
%   number of iterations to run the algorithm for.
%
%   This implements Algorithm 1 in Tan, Burrus (2019).

    if nargin < 5
        num_iters = 50;
    end
    if nargin < 4
        dc_inds = [];
    end
    
    % Convert D, b_init, and a_init into column vectors
    if size(b_init, 1) == 1
        b_init = permute(b_init, [2 1]);
    end
    if size(a_init, 1) == 1
        a_init = permute(a_init, [2 1]);
    end
    if size(D, 1) == 1
        D = permute(D, [2 1]);
    end
    nB = length(b_init) - 1;
    nA = length(a_init) - 1;

    % Prepare DFT matrices
    L = size(D,1)/2;
    % DFT matrix for b
    Wb = dftmtx(2*L);
    Wb = Wb(:,1:nB+1);
    % DFT matrix for a
    Wa = dftmtx(2*L);
    Wa = Wa(:,2:nA+1); % remove first column as well to enforce a(1)=1
    % Remove don't care regions
    Wb(dc_inds,:) = [];
    Wa(dc_inds,:) = [];
    D(dc_inds,:) = [];
    
    % Prepare initialization
    h = [b_init; a_init(2:end)];
    Bw = Wb*h(1:nB+1);
    Aw = Wa*h(nB+2:end) + ones([size(Wa,1) 1]);
    res = Bw./Aw - D;
    
    % Iterate
    best_err = inf;
    for i = 1:num_iters+1
        
        % Obtain Jacobian
        Jac = zeros(size(D,1), nB+nA+1);
        Jac(:, 1:nB+1) = Wb ./ repmat(Aw,[1 nB+1]);
        Jac(:, nB+2:end) = -(Wa.*repmat((Bw./(Aw.^2)), [1 nA]));
        % Take step
        h = real(h - (Jac \ res));
        % Flip a back into stability:
        a = real([1 h(nB+2:end)']);
        a = stabilize_poles(a);
        h(nB+2:end) = a(2:end);
        
        % Solve for freq responses
        Bw = Wb*h(1:nB+1);
        Aw = Wa*h(nB+2:end) + ones([size(Wa,1) 1]);
        % Solve for residual
        res = Bw./Aw - D;
        
        % Keep track of the best-performing filter
        new_err = norm(res)^2;
        if new_err < best_err
            best_err = new_err;
            best_h = h;
        end
    end
    b = real(best_h(1:nB+1))';
    a = real([1 best_h(nB+2:end)']);

