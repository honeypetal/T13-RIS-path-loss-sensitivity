%% Team 13 - Methodology 2 extension, Step C: additional channel parameters
% Compares SimRIS and the geometric model on parameters beyond average
% capacity (no new simulations - uses the saved channel files):
%   1. Direct-link path loss (mean, median, spread)
%   2. Rician K-factor of the RIS links (Tx->RIS and RIS->Rx)
%   3. RIS-to-direct power ratio
%   4. Outage probability P(SNR < 10 dB) with and without RIS
%   5. 5th-percentile ("cell-edge") SNR with and without RIS
%   6. Channel hardening (coefficient of variation of the channel gain)
% Run inside C:\T13_RIS (needs the SimRIS .mat files and results_M2).
clear; clc; close all;
outDir = 'results_M2X';  if ~exist(outDir, 'dir'), mkdir(outDir); end
Pt = 1; Pn = 10^((-100-30)/10);  thr_dB = 10;
cfgs  = {'InH_S1_f28', 'InH_S1_f73', 'InH_S2_f28', 'UMi_S1_f28'};
Nlist = [16 64 144 256 400];
models = {'SimRIS', 'Geometric'};

rows = cell(0, 16);
for c = 1:numel(cfgs)
    for N = Nlist
        S  = load(sprintf('%s_N%d.mat', cfgs{c}, N));
        Gm = load(fullfile('results_M2', sprintf('GEO_%s_N%d.mat', cfgs{c}, N)));
        data = {{S.H, S.G, S.D}, {Gm.Hg, Gm.Gg, Gm.Dg}};
        for m = 1:2
            [H, G, D] = deal(data{m}{:});
            st = stats(H, G, D, Pt, Pn, thr_dB);
            rows(end+1,:) = {cfgs{c}, N, models{m}, st.PL_mean, st.PL_median, st.PL_std, ...
                st.K_H_dB, st.K_G_dB, st.ratio_dB, st.out_no, st.out_opt, ...
                st.snr5_no, st.snr5_opt, st.cv_no, st.cv_opt, st.snr_med_opt}; %#ok<SAGROW>
        end
    end
    fprintf('%-11s done\n', cfgs{c});
end
T = cell2table(rows, 'VariableNames', {'Config','N','Model','DirectPL_mean_dB','DirectPL_median_dB', ...
    'DirectPL_std_dB','K_TxRIS_dB','K_RISRx_dB','RIS_to_direct_dB','Outage_noRIS','Outage_optRIS', ...
    'SNR5_noRIS_dB','SNR5_optRIS_dB','CV_noRIS','CV_optRIS','SNR_median_optRIS_dB'});
writetable(T, fullfile(outDir, 'C_channel_stats.csv'));
disp(T);

% compact summary at N = 400
S400 = T(T.N == 400, {'Config','Model','DirectPL_median_dB','DirectPL_std_dB','K_TxRIS_dB', ...
    'K_RISRx_dB','RIS_to_direct_dB','Outage_noRIS','Outage_optRIS','SNR5_noRIS_dB','SNR5_optRIS_dB'});
writetable(S400, fullfile(outDir, 'C_summary_N400.csv'));
fprintf('\n===== Channel parameters at N = 400 =====\n'); disp(S400);

col = {[0 0.45 0.74], [0.85 0.33 0.10]};  lsty = {'-o', '--s'};
% Figure 1: outage probability vs N
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:numel(cfgs)
    subplot(2,2,c); hold on; grid on;
    for m = 1:2
        t = T(strcmp(T.Config, cfgs{c}) & strcmp(T.Model, models{m}), :);
        plot(t.N, 100*t.Outage_optRIS, lsty{m}, 'Color', col{m}, 'LineWidth', 1.5);
        plot(t.N, 100*t.Outage_noRIS, ':', 'Color', col{m}, 'LineWidth', 1.5);
    end
    set(gca, 'XScale', 'log'); ylim([0 100]); xlabel('N'); ylabel(sprintf('P(SNR < %d dB)  (%%)', thr_dB));
    title(['Outage probability: ' strrep(cfgs{c}, '_', ' ')]);
    if c == 1, legend('SimRIS, opt RIS', 'SimRIS, no RIS', 'Geometric, opt RIS', 'Geometric, no RIS', 'Location', 'best', 'FontSize', 7); end
end
saveas(gcf, fullfile(outDir, 'C_outage.png'));

% Figure 2: 5th-percentile SNR vs N
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:numel(cfgs)
    subplot(2,2,c); hold on; grid on;
    for m = 1:2
        t = T(strcmp(T.Config, cfgs{c}) & strcmp(T.Model, models{m}), :);
        plot(t.N, t.SNR5_optRIS_dB, lsty{m}, 'Color', col{m}, 'LineWidth', 1.5);
        plot(t.N, t.SNR5_noRIS_dB, ':', 'Color', col{m}, 'LineWidth', 1.5);
    end
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('5th-percentile SNR (dB)');
    title(['Cell-edge SNR: ' strrep(cfgs{c}, '_', ' ')]);
    if c == 1, legend('SimRIS, opt RIS', 'SimRIS, no RIS', 'Geometric, opt RIS', 'Geometric, no RIS', 'Location', 'best', 'FontSize', 7); end
end
saveas(gcf, fullfile(outDir, 'C_snr5.png'));

% Figure 3: K-factor and RIS-to-direct ratio at N = 400
figure('Color','w', 'Position', [100 100 1000 420]);
subplot(1,2,1);
K = zeros(numel(cfgs), 4);
for c = 1:numel(cfgs)
    a = S400(strcmp(S400.Config, cfgs{c}) & strcmp(S400.Model, 'SimRIS'), :);
    b = S400(strcmp(S400.Config, cfgs{c}) & strcmp(S400.Model, 'Geometric'), :);
    K(c,:) = [a.K_TxRIS_dB, b.K_TxRIS_dB, a.K_RISRx_dB, b.K_RISRx_dB];
end
bar(categorical(strrep(cfgs,'_',' '), strrep(cfgs,'_',' ')), K); grid on;
ylabel('Median Rician K-factor (dB)');
legend('Tx-RIS SimRIS', 'Tx-RIS Geometric', 'RIS-Rx SimRIS', 'RIS-Rx Geometric', 'Location', 'best', 'FontSize', 7);
title('Rician K-factor of the RIS links (N = 400)');
subplot(1,2,2);
Rr = zeros(numel(cfgs), 2);
for c = 1:numel(cfgs)
    Rr(c,1) = S400.RIS_to_direct_dB(strcmp(S400.Config, cfgs{c}) & strcmp(S400.Model, 'SimRIS'));
    Rr(c,2) = S400.RIS_to_direct_dB(strcmp(S400.Config, cfgs{c}) & strcmp(S400.Model, 'Geometric'));
end
bar(categorical(strrep(cfgs,'_',' '), strrep(cfgs,'_',' ')), Rr); grid on;
ylabel('RIS path / direct path power (dB)'); legend('SimRIS', 'Geometric', 'Location', 'best');
title('Strength of the RIS path relative to the direct path (N = 400)');
saveas(gcf, fullfile(outDir, 'C_K_and_ratio.png'));

% Figure 4: channel hardening
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:numel(cfgs)
    subplot(2,2,c); hold on; grid on;
    for m = 1:2
        t = T(strcmp(T.Config, cfgs{c}) & strcmp(T.Model, models{m}), :);
        plot(t.N, t.CV_optRIS, lsty{m}, 'Color', col{m}, 'LineWidth', 1.5);
    end
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('CV of |h_{eff}|^2 (opt RIS)');
    title(['Channel hardening: ' strrep(cfgs{c}, '_', ' ')]);
    if c == 1, legend('SimRIS', 'Geometric', 'Location', 'best'); end
end
saveas(gcf, fullfile(outDir, 'C_hardening.png'));
fprintf('\nStep C done. Results in "%s".\n', outDir);

%% ---- local functions ----
function st = stats(H, G, D, Pt, Pn, thr_dB)
Nsym = size(H, 3); hopt = zeros(Nsym,1); hdir = hopt; ris = hopt;
KH = nan(Nsym,1); KG = KH;
for s = 1:Nsym
    hk = H(:,1,s); gk = G(1,:,s).'; dk = D(1,1,s); cc = gk .* hk;
    hopt(s) = dk + sum(cc .* exp(1i*(angle(dk) - angle(cc))));
    hdir(s) = dk;  ris(s) = sum(abs(cc));
    KH(s) = kfactor(abs(hk).^2);  KG(s) = kfactor(abs(gk).^2);
end
pl = -20*log10(abs(hdir));
snr_no  = 10*log10(Pt*abs(hdir).^2/Pn);
snr_opt = 10*log10(Pt*abs(hopt).^2/Pn);
st.PL_mean = mean(pl); st.PL_median = median(pl); st.PL_std = std(pl);
st.K_H_dB = 10*log10(median(KH, 'omitnan'));  st.K_G_dB = 10*log10(median(KG, 'omitnan'));
st.ratio_dB = 10*log10(mean(ris.^2) / mean(abs(hdir).^2));
st.out_no = mean(snr_no < thr_dB);  st.out_opt = mean(snr_opt < thr_dB);
st.snr5_no = pctl(snr_no, 5);  st.snr5_opt = pctl(snr_opt, 5);
g_no = abs(hdir).^2; g_opt = abs(hopt).^2;
st.cv_no = std(g_no)/mean(g_no);  st.cv_opt = std(g_opt)/mean(g_opt);
st.snr_med_opt = median(snr_opt);
end

function K = kfactor(p)
% Moment-based Rician K estimate from power samples across RIS elements
% (all elements share the same large-scale gain in one realization)
if numel(p) < 4, K = NaN; return; end
g = var(p) / mean(p)^2;
if g >= 1, K = 1e-3;                 % Rayleigh or worse -> K ~ 0
elseif g < 1e-4, K = 1e4;            % essentially pure LOS (cap at 40 dB)
else, K = sqrt(1-g) / (1 - sqrt(1-g));
end
end

function v = pctl(x, p)
% percentile without the Statistics Toolbox
x = sort(x(:)); n = numel(x);
v = x(max(1, min(n, round(p/100*n))));
end
