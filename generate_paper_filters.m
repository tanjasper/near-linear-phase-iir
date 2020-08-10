clear all; clc;
addpath('fns');

% Parameters
L = 1024; % number of freq. samples (L+1 samples from 0 to pi)
fp1 = 0.125; % passband edge from 0 to 0.5
fs1 = fp1 + 0.01; % stopband edge from 0 to 0.5
num_iters = 3000;
fir_zeros = 20;

% Generate ideal magnitude response
[D_ideal, dc_inds] = generate_ideal_mag_response([0 fp1*2 fs1*2 1], [1 1 0 0], L, 'whole');
pb_edge = dc_inds(1)-1; % for LPF, the passband edge is right before the start of the don't care region
% Generate a response with FIR phase
b_fir = firls(fir_zeros, [0 fp1*2 fs1*2 1], [1 1 0 0]);
D_fir = abs(D_ideal) .* exp(1i*angle(freqz(b_fir, 1, L*2, 'whole')));

save_name = sprintf('data/filters_diff_zeros_poles_L%d_fp0p%03d_fs0p%03d_iters%d.mat', L, fp1*1000, fs1*1000, num_iters);

zero_vals = [2 4 6 8 10 12 14 16 18 20];
poles = [2 4 6 8 10 12 14 16 18 20];

l2_errs = zeros([max(zero_vals)*2 max(zero_vals) max(poles)]);
mean_l2_errs = zeros([max(zero_vals)*2 max(zero_vals) max(poles)]);
phase_r2s = zeros([max(zero_vals)*2 max(zero_vals) max(poles)]);
grp_delays = zeros([max(zero_vals)*2 max(zero_vals) max(poles)]);

% 5, 7, 9
for p = 1:length(poles)  % for # of poles
    nA = poles(p);
    tic
   for zz = 1:length(zero_vals)  % for # of zeros
        tic
        nB = zero_vals(zz);
        b_est = zeros([nB+1 max(zero_vals)*2]);
        a_est = zeros([nA+1 max(zero_vals)*2]);
        parfor z = 1:max(zero_vals)*2  % for diff desired responses
            fprintf('p: %d out of %d, zz: %d out of %d, z: %d out of %d\n', p, length(poles), zz, length(zero_vals), z, max(zero_vals)*2);
            b_firs{z} = firls(z, [0 fp1*2 fs1*2 1], [1 1 0 0]);
            D_fir = abs(D_ideal) .* exp(1i*angle(freqz(b_firs{z}, 1, L*2, 'whole')));

            % perform optimization
            [b_init, a_init] = smb_mode1(D_fir, nB, nA);
            [b_temp, a_temp] = gauss_newton_iir(D_fir,b_init,a_init,dc_inds,num_iters);

            % parse optimization
            b_est(:,z) = b_temp';
            a_est(:,z) = a_temp';

            % Check errors
            [l2_errs(z,nB,nA), mean_l2_errs(z, nB, nA)] = mag_err(b_est(:,z), a_est(:,z), D_fir(1:L+1), dc_inds);
            [phase_r2s(z,nB,nA), grp_delays(z,nB,nA)] = lin_phase_r2(b_est(:,z), a_est(:,z), L, 1:pb_edge, false);
        end
        b_all{nB,nA} = b_est; a_all{nB,nA} = a_est;
        save(save_name, 'b_all', 'a_all', 'l2_errs', 'phase_r2s', 'grp_delays');
    end
    save('check2.mat', 'b_all', 'a_all', 'l2_errs', 'phase_r2s', 'grp_delays');
    toc
end