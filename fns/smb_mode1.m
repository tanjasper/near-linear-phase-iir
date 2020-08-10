function [ b, a ] = smb_mode1( D, nB, nA, X )
% [ b, a ] = smb_mode1( D, nB, nA )
%   Performs Steiglitz-McBride mode 1 IIR design
%   The method is simply a least squares minimization in time domain
%   D -- frequency response on uniform frequency samples around whole unit
%   circle
%   nB -- # of poles
%   nA -- # of zeros
%   X -- optional input frequency response (with the objective that the
%   input filtered by the target filter yields the output D)

num_samples = 1024;

% Step 1: obtain time domain of desired response
y = ifft(D, num_samples);
if nargin < 4
    x = [1];
else
    x = ifft(X, num_samples);
end

% Step 2: Obtain q vectors via convolution matrices
Xconv = convmtx(x, num_samples*2);
Yconv = convmtx([0; y], num_samples*2);
Xconv = Xconv(1:num_samples*2,1:nB+1);
Yconv = Yconv(1:num_samples*2,1:nA);
q_vecs = [Xconv Yconv]; % each row of this matrix is q_j'

% Step 3: Obtain Q and c
Q = zeros(nB+nA+1);
c = zeros([nB+nA+1 1]);
for i = 1:num_samples
    Q = Q + q_vecs(i,:)' * q_vecs(i,:);
    c = c + y(i) * q_vecs(i,:)';
end

% Step 4: Solve inverse problem
delta = Q \ c;
b = real(delta(1:nB+1))';
a = real([1; -delta(nB+2:end)])';
a = stabilize_poles(a);

end

