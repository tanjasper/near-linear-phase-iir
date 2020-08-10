%   Generates all the paper's results
%   Can take order of 10s. to run
%   Generates Fig. 1 and Fig. 2 from paper and displays Table 1 results

%% Load generated filters (required for all subsequent sections)

clear all; clc;
addpath('fns');

% Load generated IIR filters (result of generate_paper_filters.m)
load('data/filters_diff_zeros_poles_L1024_fp0p125_fs0p135_iters3000.mat');
L = 1024;  % number of freq. samples (L+1 samples from 0 to pi)
fp1 = 0.125;  % passband edge freq (in units of 2pi)
fs1 = 0.135;  % stopband edge freq (in units of 2pi)

% Generate ideal magnitude response
[D_ideal, dc_inds] = generate_ideal_mag_response([0 fp1*2 fs1*2 1], [1 1 0 0], L, 'whole');

% The number of zeros for which we have filters
zero_vals = [2 4 6 8 10 12 14 16 18 20];

% Average out l2 errors from the loaded mat file
num_mag_freqs = L+1 - sum(dc_inds <= (L));
mean_l2_errs = l2_errs / num_mag_freqs;


%% Figure 1: Errors & phase R2 values for diff poles, zeros, target delays

mean_fir_err = zeros(40, 1);
for i = 1:40
    [~, mean_fir_err(i)] = mag_err(firls(i, [0 fp1*2 fs1*2 1], [1 1 0 0])', 1, D_ideal(1:L+1), dc_inds);
end

figure
for i = 1:9
    subplot(2,5,i)
    imagesc(log(l2_errs(:,2:2:end,2*i)/(L+1)), [-8.5 -3.5]), colormap gray
    ylabel('Target group delay');
    xlabel('# zeros');
    title(sprintf('%d poles', i*2))
    set(gca, 'XTick', 2:2:10, 'XTickLabel', 4:4:20);
    set(gca, 'YTick', 10:10:40, 'YTickLabel', 5:5:20);
    colormap gray
    cbh = colorbar;
    % for log
    cbh.Ticks = [log(0.0003), log(0.001), log(0.005), log(0.03)];
    cbh.TickLabels = [3e-4, 1e-3, 5e-3, 3e-2];
end
subplot(2,5,10)
imagesc(log(mean_fir_err), [-8.5 -3.5]), colormap gray
ylabel('Target group delay');
title('FIR filters')
set(gca, 'YTick', 10:10:40, 'YTickLabel', 5:5:20);
set(gca,'xtick',[])
cbh = colorbar;
% for log
cbh.Ticks = [log(0.0003), log(0.001), log(0.005), log(0.03)];
cbh.TickLabels = [3e-4, 1e-3, 5e-3, 3e-2];

figure
for i = 1:9
    subplot(2,5,i)
    imagesc(phase_r2s(:,2:2:end,2*i), [0.98 1]), colormap gray
    ylabel('Target group delay');
    xlabel('# zeros');
    title(sprintf('%d poles', i*2))
    set(gca, 'XTick', 2:2:10, 'XTickLabel', 4:4:20);
    set(gca, 'YTick', 10:10:40, 'YTickLabel', 5:5:20);
    colormap gray
    colorbar
end


%% Table 1: Find min zeros to achieve l2 error

fir_zeros = 80;

[~, mean_fir_err] = mag_err(firls(fir_zeros, [0 fp1*2 fs1*2 1], [1 1 0 0])', 1, D_ideal(1:L+1), dc_inds);

num_poles_all = [2 4 6 8 10 12];
mean_l2_val = mean_fir_err;

fprintf('\n\n========== FIR error: %.6f ==========\n\n', mean_fir_err);
for p = 1:length(num_poles_all)
    num_poles = num_poles_all(p);
    curr_l2_errs = mean_l2_errs(:,:,num_poles);
    curr_phase_r2s = phase_r2s(:,:,num_poles);
    curr_grp_delays = grp_delays(:,:,num_poles);

    for zz = 1:length(zero_vals)
        temp_l2_errs = mean_l2_errs(:,zero_vals(zz), num_poles);
        temp_phase_r2s = phase_r2s(:,zero_vals(zz), num_poles);
        temp_grp_delays = grp_delays(:,zero_vals(zz), num_poles);
        if min(temp_l2_errs) > mean_l2_val % does not achieve the desired l2 value
            continue
        else
            % smallest l2 error
            [~,idx_l2] = min(temp_l2_errs);
            num_zeros = zero_vals(zz);
            fprintf('========== Results for %d poles ==========\n', num_poles);
            fprintf('Min l2err -- initialization zeros: %d, l2err: %.6f, phase r2: %.6f, grp_delay: %.3f\n', idx_l2, temp_l2_errs(idx_l2), temp_phase_r2s(idx_l2), temp_grp_delays(idx_l2));
            fprintf('# zeros: %d, Delays: %d, Unique coefficients: %d, Multiplies: %d, Adds: %d\n', zero_vals(zz), max(num_zeros, num_poles), ceil((num_zeros+1)/2)+num_poles, num_zeros+num_poles+1, num_zeros+num_poles);
            coeffs_l2(p) = ceil((num_zeros+1)/2)+num_poles;
            mults_l2(p) = num_zeros+num_poles+1;
            % best phase value
            [~,idx_r2s] = sort(temp_phase_r2s,'descend');
            for aa = 1:length(idx_r2s)
                if temp_l2_errs(idx_r2s(aa)) <= mean_l2_val
                    num_zeros = zero_vals(zz);
                    idx_phase = idx_r2s(aa);
                    fprintf('Most linear phase -- initialization zeros: %d, l2err: %.6f, phase r2: %.6f, grp_delay: %.3f\n', idx_phase, temp_l2_errs(idx_phase), temp_phase_r2s(idx_phase), temp_grp_delays(idx_phase));
                    fprintf('# zeros: %d, Delays: %d, Unique coefficients: %d, Multiplies: %d, Adds: %d\n', num_zeros, max(num_zeros, num_poles), ceil((num_zeros+1)/2)+num_poles, num_zeros+num_poles+1, num_zeros+num_poles);
                    coeffs_phase(p) = ceil((num_zeros+1)/2)+num_poles;
                    mults_phase(p) = num_zeros+num_poles+1;
                    break
                end
            end
            break
        end
    end
end


%% Figure 2: Plot with multiple FIR filters and a better IIR filter for each

% helper variables
num_tried_zeros = size(l2_errs, 2);
num_tried_poles = size(l2_errs, 3);

% ignore entries without results
l2_errs(l2_errs == 0) = inf;

% Parameters
phase_r2_thresh = 0.99;

% Which IIR configurations pass the phase threshold?
qualifying_entries_phase = find(phase_r2s >= phase_r2_thresh);
unqualifying_entries_phase = find(phase_r2s < phase_r2_thresh);

% For desired magnitude response, we will use FIR filters performance
[D_ideal, dc_inds] = generate_ideal_mag_response([0 fp1*2 fs1*2 1], [1 1 0 0], L, 'whole');
fir_zeros = 8:2:124;

% how many multiplies and coefficients for each IIR configuration?
num_mults = zeros(size(l2_errs, 2), size(l2_errs, 3));
num_coeffs = zeros(size(l2_errs, 2), size(l2_errs, 3));
for nB = 1:size(l2_errs, 2)
    for nA = 1:size(l2_errs, 3)
        num_mults(nB, nA) = nB + nA + 1;
        num_coeffs(nB, nA) = ceil((nB+1)/2) + nA;
    end
end

% For each FIR error, find the possible phase r2s achievable
best_phase_r2s = zeros(size(phase_r2s, 2), size(phase_r2s, 3), length(fir_zeros));
for i = 1:length(fir_zeros)
    
    % FIR error
    [~, mean_fir_err(i)] = mag_err(firls(fir_zeros(i), [0 fp1*2 fs1*2 1], [1 1 0 0])', 1, D_ideal(1:L+1), dc_inds);
    
    % find IIR configurations that achieves this error or lower
    qualifying_entries_mag = find(mean_l2_errs <= mean_fir_err(i));
    qualifying_entries = intersect(qualifying_entries_mag, qualifying_entries_phase); % configurations that satisfy both magnitude and phase requirements
    unqualifying_entries_mag = find(mean_l2_errs > mean_fir_err(i));
    unqualifying_entries = union(unqualifying_entries_mag, unqualifying_entries_phase);
    
    % ignore the phase initializations
    unqualifying_flags = zeros(size(mean_l2_errs));
    unqualifying_flags(unqualifying_entries) = 1;
    unqualifying_flags = squeeze(min(unqualifying_flags));
    % empty entries are unqualifying
    unqualifying_flags(1:2:end, :) = 1;
    unqualifying_flags(:, 1:2:end) = 1;
    
    % which # coeffs satisfy magnitude and phase requirements?
    curr_num_coeffs = num_coeffs;
    curr_num_coeffs(unqualifying_flags==1) = inf;  % did not qualify, set # coefficients to inf to ignore
    best_coeffs(i) = min(curr_num_coeffs(:));  % minimum number of coefficients
    best_inds = find(curr_num_coeffs == best_coeffs(i));  % indices of configurations that achieved the minimum # coeffs
    temp_achieving_zeros = 0;
    temp_achieving_poles = 0;
    for a = 1:length(best_inds)
        curr_ind = best_inds(a);
        temp_achieving_zeros(a) = mod(curr_ind-1, num_tried_zeros) + 1;
        temp_achieving_poles(a) = ceil(curr_ind / num_tried_zeros);
    end
    achieving_zeros_coeffs{i} = temp_achieving_zeros;
    achieving_poles_coeffs{i} = temp_achieving_poles;
    
    % do exact same thing for # mults
    curr_num_mults = num_mults;
    curr_num_mults(unqualifying_flags==1) = inf;
    best_mults(i) = min(curr_num_mults(:));  % minimum number of coefficients
    best_inds = find(curr_num_mults == best_mults(i));  % indices of configurations that achieved the minimum # coeffs
    temp_achieving_zeros = 0;
    temp_achieving_poles = 0;
    temp_achieving_mags = 0;
    for a = 1:length(best_inds)
        curr_ind = best_inds(a);
        temp_achieving_zeros(a) = mod(curr_ind-1, num_tried_zeros) + 1;
        temp_achieving_poles(a) = ceil(curr_ind / num_tried_zeros);
    end
    achieving_zeros_mults{i} = temp_achieving_zeros;
    achieving_poles_mults{i} = temp_achieving_poles;
    
end

% What are the magnitudes of those qualifying entries?
for i = 1:length(fir_zeros)
    min_l2_err = inf;
    for a = 1:length(achieving_zeros_mults{i})
        nz = achieving_zeros_mults{i}(a);
        np = achieving_poles_mults{i}(a);
        curr_l2_errs = mean_l2_errs(:, nz, np);
        curr_phases = phase_r2s(:, nz, np);
        curr_l2_errs(curr_l2_errs == 0) = inf;
        curr_l2_errs(curr_phases < phase_r2_thresh) = inf;
        curr_l2_err = min(curr_l2_errs);
        if curr_l2_err < min_l2_err
            min_l2_err = curr_l2_err;
            best_zeros_mults(i) = nz;
            best_poles_mults(i) = np;
        end
    end
    best_l2_errs(i) = min_l2_err;
end

figure
plot(fir_zeros, mean_fir_err, 'linestyle', '--', 'linewidth', 2), hold on
plot(fir_zeros, best_l2_errs, 'linestyle', '-', 'linewidth', 2), hold on
yl1 = ylabel('Magnitude mean L_2 error');
l1 = legend('FIR', 'IIR ($R^2 \geq 0.99$)');
xl1 = xlabel('FIR order');
xlim([8 124])
t1 = title('Magnitude Error');
set(gca, 'FontSize', 16)
figure
plot(fir_zeros, fir_zeros+1, 'linestyle', '--', 'linewidth', 2), hold on
plot(fir_zeros, best_mults, 'linestyle', '-', 'linewidth', 2)
yl2 = ylabel('Number of multiplies');
xl2 = xlabel('FIR order');
l2 = legend('FIR', 'IIR ($R^2 \geq 0.99$)');
xlim([8 124])
t2 = title('Number of Multiplies');
set(gca, 'FontSize', 16)

% convert all text to latex
set(yl1, 'string', 'Magnitude mean $\ell_2$ error', 'interpreter', 'latex');
set(yl2, 'string', 'Number of multiplies', 'interpreter', 'latex');
set(xl1, 'string', 'FIR order', 'interpreter', 'latex');
set(xl2, 'string', 'FIR order', 'interpreter', 'latex');
set(l1, 'Interpreter', 'latex');
set(l2, 'Interpreter', 'latex');
set(t1, 'string', '\textbf{Magnitude Error}', 'interpreter', 'latex');
set(t2, 'string', '\textbf{Number of Multiplies}', 'interpreter', 'latex');
set(yl1, 'FontSize', 18)
set(yl2, 'FontSize', 18)
set(xl1, 'FontSize', 18)
set(xl2, 'FontSize', 18)
set(l1, 'FontSize', 16)
set(l2, 'FontSize', 16)
set(t1, 'FontSize', 18)
set(t2, 'FontSize', 18)

