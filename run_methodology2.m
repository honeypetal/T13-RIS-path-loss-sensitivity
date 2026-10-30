%% Team 13 - Methodology 2: Cross-simulator validation (SimRIS vs geometric model)
% Put this file and geo_ris_channel.m in C:\T13_RIS (the same folder as
% your 20 SimRIS .mat files from Methodology 1), then press Run.
% Nothing needs to be clicked in SimRIS - the second model runs by itself.
clear; clc; close all;

Pt_dBm    = 30;      % same as Methodology 1
noise_dBm = -100;
Nsym      = 100;     % same number of realizations as SimRIS
Nlist     = [16 64 144 256 400];
outDir    = 'results_M2';
if ~exist(outDir, 'dir'), mkdir(outDir); end
rng(2026);           % makes the geometric model repeatable

% name, environment, scenario, frequency, Tx, Rx, RIS  (same as SimRIS runs)
cfg = {
 'InH_S1_f28', 'InH', 1, 28, [0 25 2],  [38 48 1], [40 50 2];
 'InH_S1_f73', 'InH', 1, 73, [0 25 2],  [38 48 1], [40 50 2];
 'InH_S2_f28', 'InH', 2, 28, [0 25 2],  [70 35 1], [70 30 2];
 'UMi_S1_f28', 'UMi', 1, 28, [0 25 20], [80 75 1], [70 85 10] };

Pt = 10^((Pt_dBm-30)/10);  Pn = 10^((noise_dBm-30)/10);

%% 1) Run both models for every configuration and N
rows = cell(0, 12);
snrStore = struct();
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for N = Nlist
        simFile = sprintf('%s_N%d.mat', name, N);
        if ~exist(simFile, 'file')
            error('Missing SimRIS file %s - run this script inside C:\\T13_RIS', simFile);
        end
        S = load(simFile);
        [Cs, Ps] = evaluate(S.H, S.G, S.D, Pt, Pn);

        [Hg, Gg, Dg] = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, Nsym);
        [Cg, Pg] = evaluate(Hg, Gg, Dg, Pt, Pn);

        rows(end+1,:) = {name, N, Cs.noRIS, Cg.noRIS, Cs.opt, Cg.opt, ...
            Cs.opt - Cs.noRIS, Cg.opt - Cg.noRIS, Ps.noRIS, Pg.noRIS, Ps.opt, Pg.opt}; %#ok<SAGROW>
        fprintf('%-12s N=%3d | no RIS: SimRIS %5.2f  Geo %5.2f | opt RIS: SimRIS %5.2f  Geo %5.2f bps/Hz\n', ...
            name, N, Cs.noRIS, Cg.noRIS, Cs.opt, Cg.opt);
        if N == 400
            snrStore.(name) = struct('sim_opt', Cs.snr_opt, 'sim_no', Cs.snr_no, ...
                                     'geo_opt', Cg.snr_opt, 'geo_no', Cg.snr_no);
        end
        save(fullfile(outDir, sprintf('GEO_%s_N%d.mat', name, N)), 'Hg', 'Gg', 'Dg');
    end
end
T = cell2table(rows, 'VariableNames', {'Config','N', ...
    'C_noRIS_SimRIS','C_noRIS_Geo','C_opt_SimRIS','C_opt_Geo', ...
    'Gain_SimRIS','Gain_Geo','Prx_noRIS_SimRIS_dBm','Prx_noRIS_Geo_dBm', ...
    'Prx_opt_SimRIS_dBm','Prx_opt_Geo_dBm'});
writetable(T, fullfile(outDir, 'M2_comparison.csv'));
disp(T);

%% 2) Agreement metrics (SimRIS taken as reference)
metric = @agreement;
qty = {'C_opt', 'C_opt_SimRIS', 'C_opt_Geo';
       'C_noRIS', 'C_noRIS_SimRIS', 'C_noRIS_Geo';
       'Gain', 'Gain_SimRIS', 'Gain_Geo';
       'Prx_opt_dBm', 'Prx_opt_SimRIS_dBm', 'Prx_opt_Geo_dBm'};
groups = [cfg(:,1).', {'All_indoor', 'All'}];
mrows = cell(0, 6);
for g = 1:numel(groups)
    switch groups{g}
        case 'All',        sel = true(height(T),1);
        case 'All_indoor', sel = startsWith(T.Config, 'InH');
        otherwise,         sel = strcmp(T.Config, groups{g});
    end
    for q = 1:size(qty,1)
        v = metric(T.(qty{q,2})(sel), T.(qty{q,3})(sel));
        mrows(end+1,:) = {groups{g}, qty{q,1}, v(1), v(2), v(3), v(4)}; %#ok<SAGROW>
    end
end
Mt = cell2table(mrows, 'VariableNames', {'Group','Quantity','RMSE','MAE','R2','Pearson_r'});
writetable(Mt, fullfile(outDir, 'M2_metrics.csv'));
fprintf('\n===== Agreement metrics (SimRIS = reference) =====\n');
disp(Mt);

%% 3) Plot: capacity vs N, both simulators, one figure per configuration
for c = 1:size(cfg,1)
    t = T(strcmp(T.Config, cfg{c,1}), :);
    figure('Name', cfg{c,1}, 'Color', 'w'); hold on; grid on;
    plot(t.N, t.C_opt_SimRIS,   '-o',  'LineWidth', 1.6, 'Color', [0 0.45 0.74]);
    plot(t.N, t.C_opt_Geo,      '--o', 'LineWidth', 1.6, 'Color', [0.30 0.75 0.93]);
    plot(t.N, t.C_noRIS_SimRIS, '-s',  'LineWidth', 1.6, 'Color', [0.85 0.33 0.10]);
    plot(t.N, t.C_noRIS_Geo,    '--s', 'LineWidth', 1.6, 'Color', [0.93 0.69 0.13]);
    set(gca, 'XScale', 'log'); xlabel('Number of RIS elements N'); ylabel('Ergodic capacity (bps/Hz)');
    title(['SimRIS vs geometric model: ' strrep(cfg{c,1}, '_', ' ')]);
    legend('Optimised RIS - SimRIS', 'Optimised RIS - Geometric', 'No RIS - SimRIS', 'No RIS - Geometric', 'Location', 'northwest');
    saveas(gcf, fullfile(outDir, ['M2_Cap_vs_N_' cfg{c,1} '.png']));
end

%% 4) Plot: scatter agreement (Geo vs SimRIS) with y = x line
figure('Name', 'Scatter', 'Color', 'w');
pairs = {'C_opt_SimRIS','C_opt_Geo','Optimised RIS capacity'; 'C_noRIS_SimRIS','C_noRIS_Geo','No-RIS capacity'};
for p = 1:2
    subplot(1,2,p); hold on; grid on; axis square;
    mk = {'o','s','^','d'};
    for c = 1:size(cfg,1)
        t = T(strcmp(T.Config, cfg{c,1}), :);
        plot(t.(pairs{p,1}), t.(pairs{p,2}), mk{c}, 'MarkerSize', 7, 'LineWidth', 1.4, ...
             'DisplayName', strrep(cfg{c,1}, '_', ' '));
    end
    lims = [0, ceil(max([T.(pairs{p,1}); T.(pairs{p,2})]))+1];
    plot(lims, lims, 'k--', 'DisplayName', 'y = x'); xlim(lims); ylim(lims);
    v = metric(T.(pairs{p,1}), T.(pairs{p,2}));
    xlabel('SimRIS (bps/Hz)'); ylabel('Geometric model (bps/Hz)');
    title(sprintf('%s\nRMSE = %.2f, R^2 = %.2f', pairs{p,3}, v(1), v(3)));
    legend('Location', 'southeast', 'FontSize', 7);
end
set(gcf, 'Position', [100 100 1000 450]);
saveas(gcf, fullfile(outDir, 'M2_scatter.png'));

%% 5) Plot: RIS capacity gain, both simulators
figure('Name', 'Gain', 'Color', 'w'); hold on; grid on;
col = lines(size(cfg,1));
for c = 1:size(cfg,1)
    t = T(strcmp(T.Config, cfg{c,1}), :);
    plot(t.N, t.Gain_SimRIS, '-o',  'Color', col(c,:), 'LineWidth', 1.6, 'DisplayName', [strrep(cfg{c,1},'_',' ') ' SimRIS']);
    plot(t.N, t.Gain_Geo,    '--s', 'Color', col(c,:), 'LineWidth', 1.6, 'DisplayName', [strrep(cfg{c,1},'_',' ') ' Geometric']);
end
set(gca, 'XScale', 'log'); xlabel('Number of RIS elements N'); ylabel('Capacity gain over no-RIS (bps/Hz)');
title('RIS capacity gain: SimRIS (solid) vs geometric model (dashed)');
legend('Location', 'northwest', 'FontSize', 7);
saveas(gcf, fullfile(outDir, 'M2_gain_compare.png'));

%% 6) Plot: RMSE per configuration
figure('Name', 'RMSE', 'Color', 'w');
R = zeros(size(cfg,1), 3);
for c = 1:size(cfg,1)
    for q = 1:3
        R(c,q) = Mt.RMSE(strcmp(Mt.Group, cfg{c,1}) & strcmp(Mt.Quantity, qty{q,1}));
    end
end
bar(categorical(strrep(cfg(:,1), '_', ' ')), R); grid on;
ylabel('RMSE (bps/Hz)'); legend('Optimised RIS capacity', 'No-RIS capacity', 'Capacity gain', 'Location', 'northwest');
title('Cross-simulator disagreement (RMSE) per configuration');
saveas(gcf, fullfile(outDir, 'M2_rmse_bars.png'));

%% 7) Plot: SNR CDFs at N = 400 (distribution-level comparison)
figure('Name', 'CDF', 'Color', 'w');
for c = 1:size(cfg,1)
    s = snrStore.(cfg{c,1});
    subplot(2,2,c); hold on; grid on;
    ecdfplot(s.sim_opt, '-',  [0 0.45 0.74]);   ecdfplot(s.geo_opt, '--', [0.30 0.75 0.93]);
    ecdfplot(s.sim_no,  '-',  [0.85 0.33 0.10]); ecdfplot(s.geo_no,  '--', [0.93 0.69 0.13]);
    xlabel('SNR (dB)'); ylabel('CDF'); title([strrep(cfg{c,1}, '_', ' ') ', N = 400']);
    if c == 1
        legend('Opt RIS - SimRIS', 'Opt RIS - Geo', 'No RIS - SimRIS', 'No RIS - Geo', 'Location', 'southeast', 'FontSize', 6);
    end
end
set(gcf, 'Position', [100 100 900 650]);
saveas(gcf, fullfile(outDir, 'M2_SNR_CDF_N400.png'));

fprintf('\nDone. Tables and figures are in the "%s" folder.\n', outDir);

%% ---------------- local functions ----------------
function [C, P] = evaluate(H, G, D, Pt, Pn)
% Capacity (bps/Hz), SNR (dB) and received power (dBm) for one dataset
Nsym = size(H, 3);
hopt = zeros(Nsym,1); hrnd = hopt; hdir = hopt;
for s = 1:Nsym
    hk = H(:,1,s); gk = G(1,:,s).'; dk = D(1,1,s);
    cc = gk .* hk;
    hopt(s) = dk + sum(cc .* exp(1i*(angle(dk) - angle(cc))));
    hrnd(s) = dk + sum(cc .* exp(1i*2*pi*rand(numel(cc),1)));
    hdir(s) = dk;
end
snr = @(h) Pt*abs(h).^2/Pn;
C.opt   = mean(log2(1 + snr(hopt)));
C.rand  = mean(log2(1 + snr(hrnd)));
C.noRIS = mean(log2(1 + snr(hdir)));
C.snr_opt = 10*log10(snr(hopt));
C.snr_no  = 10*log10(snr(hdir));
P.opt   = 10*log10(mean(Pt*abs(hopt).^2)) + 30;
P.noRIS = 10*log10(mean(Pt*abs(hdir).^2)) + 30;
end

function v = agreement(ref, x)
% [RMSE, MAE, R^2, Pearson r] of x against the reference ref
ref = ref(:); x = x(:);
r   = corrcoef(ref, x);
v   = [sqrt(mean((x-ref).^2)), mean(abs(x-ref)), ...
       1 - sum((x-ref).^2)/sum((ref-mean(ref)).^2), r(1,2)];
end

function ecdfplot(x, ls, colr)
x = sort(x(:)); y = (1:numel(x)).' / numel(x);
stairs(x, y, ls, 'LineWidth', 1.5, 'Color', colr);
end
