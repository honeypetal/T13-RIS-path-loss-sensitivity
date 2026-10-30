%% Team 13 - Methodology 2 extension, Step A: ablation study
% Which modelling assumption causes the SimRIS-vs-geometric disagreement?
% The geometric model is made more SimRIS-like ONE assumption at a time:
%   V0  baseline geometric model (3GPP TR 38.901 everywhere)
%   V1  + SimRIS's path-loss model (5GCM close-in, SimRIS's n and sigma, all links)
%   V2  + SimRIS's LOS-probability model (5GCM, 3-D distance)
% Whatever gap is left after V2 is due to SimRIS's scattering-cluster
% multipath model, which the geometric model deliberately does not have.
% Run inside C:\T13_RIS (needs the 20 SimRIS .mat files).
clear; clc; close all;
outDir = 'results_M2X';  if ~exist(outDir, 'dir'), mkdir(outDir); end
Pt = 1; Pn = 10^((-100-30)/10);  Nsym = 100;  Nlist = [16 64 144 256 400];

cfg = { 'InH_S1_f28', 'InH', 1, 28, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S1_f73', 'InH', 1, 73, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S2_f28', 'InH', 2, 28, [0 25 2],  [70 35 1], [70 30 2];
        'UMi_S1_f28', 'UMi', 1, 28, [0 25 20], [80 75 1], [70 85 10] };

% SimRIS's own large-scale parameters (from SimRIS_GUI.m, 5GCM model).
% NLOS exponent includes SimRIS's frequency term n*(1 + b*(f - 24.2)/24.2), b = 0.06 (InH).
variants = {'V0 baseline (3GPP)', 'V1 + SimRIS path loss', 'V2 + SimRIS LOS probability'};
rows = cell(0, 9);
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for N = Nlist
        S = load(sprintf('%s_N%d.mat', name, N));
        ref = evaluate(S.H, S.G, S.D, Pt, Pn);
        for v = 1:3
            rng(100*c + N + v);                       % repeatable
            switch v
                case 1, cal = [];
                case 2, cal = simris_pl(env, f);
                case 3, cal = simris_pl(env, f); cal.plos = '5gcm';
            end
            [H, G, D] = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, Nsym, cal);
            r = evaluate(H, G, D, Pt, Pn);
            rows(end+1,:) = {name, N, variants{v}, ref.noRIS, r.noRIS, ref.opt, r.opt, ...
                             ref.opt-ref.noRIS, r.opt-r.noRIS}; %#ok<SAGROW>
        end
        fprintf('%-11s N=%3d done\n', name, N);
    end
end
T = cell2table(rows, 'VariableNames', {'Config','N','Variant','C_noRIS_SimRIS','C_noRIS_Geo', ...
    'C_opt_SimRIS','C_opt_Geo','Gain_SimRIS','Gain_Geo'});
writetable(T, fullfile(outDir, 'A_ablation_all.csv'));

% RMSE vs SimRIS for every variant, per configuration and overall
groups = [cfg(:,1).', {'All'}];
qs = {'C_noRIS', 'C_opt', 'Gain'};
mrows = cell(0, 6);
for g = 1:numel(groups)
    for v = 1:3
        sel = strcmp(T.Variant, variants{v});
        if ~strcmp(groups{g}, 'All'), sel = sel & strcmp(T.Config, groups{g}); end
        e = zeros(1,3);
        for q = 1:3
            e(q) = sqrt(mean((T.([qs{q} '_Geo'])(sel) - T.([qs{q} '_SimRIS'])(sel)).^2));
        end
        mrows(end+1,:) = {groups{g}, variants{v}, e(1), e(2), e(3), NaN}; %#ok<SAGROW>
    end
end
M = cell2table(mrows, 'VariableNames', {'Group','Variant','RMSE_C_noRIS','RMSE_C_opt','RMSE_Gain','Gap_explained_pct'});
for g = 1:numel(groups)
    i0 = find(strcmp(M.Group, groups{g}) & strcmp(M.Variant, variants{1}));
    for v = 1:3
        i = find(strcmp(M.Group, groups{g}) & strcmp(M.Variant, variants{v}));
        M.Gap_explained_pct(i) = 100*(M.RMSE_C_noRIS(i0) - M.RMSE_C_noRIS(i))/M.RMSE_C_noRIS(i0);
    end
end
writetable(M, fullfile(outDir, 'A_ablation_summary.csv'));
fprintf('\n===== Ablation: RMSE against SimRIS (bps/Hz) =====\n'); disp(M);

% Figure 1: RMSE of no-RIS capacity per variant and configuration
figure('Color','w', 'Position', [100 100 900 420]);
B = zeros(numel(groups), 3);
for g = 1:numel(groups)
    for v = 1:3
        B(g,v) = M.RMSE_C_noRIS(strcmp(M.Group, groups{g}) & strcmp(M.Variant, variants{v}));
    end
end
bar(categorical(strrep(groups,'_',' '), strrep(groups,'_',' ')), B); grid on;
ylabel('RMSE vs SimRIS, no-RIS capacity (bps/Hz)'); legend(variants, 'Location', 'northwest');
title('Ablation: how much of the direct-link disagreement each assumption explains');
saveas(gcf, fullfile(outDir, 'A_ablation_noRIS_rmse.png'));

% Figure 2: same for optimised-RIS capacity and gain (overall)
figure('Color','w');
i = strcmp(M.Group, 'All');
bar(categorical({'No-RIS capacity','Optimised RIS capacity','RIS gain'}, {'No-RIS capacity','Optimised RIS capacity','RIS gain'}), ...
    [M.RMSE_C_noRIS(i) M.RMSE_C_opt(i) M.RMSE_Gain(i)]'); grid on;
ylabel('RMSE vs SimRIS (bps/Hz), all 20 points'); legend(variants, 'Location', 'northeast');
title('Ablation summary (all configurations)');
saveas(gcf, fullfile(outDir, 'A_ablation_summary.png'));

% Figure 3: no-RIS capacity vs N, SimRIS and the three variants
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:size(cfg,1)
    subplot(2,2,c); hold on; grid on;
    t = T(strcmp(T.Config, cfg{c,1}) & strcmp(T.Variant, variants{1}), :);
    plot(t.N, t.C_noRIS_SimRIS, 'k-o', 'LineWidth', 2);
    mk = {'--s', '--^', '--d'};
    for v = 1:3
        t = T(strcmp(T.Config, cfg{c,1}) & strcmp(T.Variant, variants{v}), :);
        plot(t.N, t.C_noRIS_Geo, mk{v}, 'LineWidth', 1.4);
    end
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('No-RIS capacity (bps/Hz)');
    title(strrep(cfg{c,1}, '_', ' '));
    if c == 1, legend([{'SimRIS'}, variants], 'Location', 'best', 'FontSize', 7); end
end
saveas(gcf, fullfile(outDir, 'A_ablation_noRIS_vs_N.png'));
fprintf('\nStep A done. Results in "%s".\n', outDir);

%% ---- local functions ----
function p = simris_pl(env, f)
% SimRIS's own large-scale parameters (5GCM close-in model, from SimRIS_GUI.m).
% InH NLOS exponent includes SimRIS's frequency term n*(1 + 0.06*(f - 24.2)/24.2).
if strcmp(env, 'InH')
    p = struct('n_LOS', 1.73, 'sigma_LOS', 3.02, ...
               'n_NLOS', 3.19*(1 + 0.06*(f - 24.2)/24.2), 'sigma_NLOS', 8.29);
else
    p = struct('n_LOS', 1.98, 'sigma_LOS', 3.10, 'n_NLOS', 3.19, 'sigma_NLOS', 8.20);
end
p.apply_ris = true;
end

function C = evaluate(H, G, D, Pt, Pn)
Nsym = size(H, 3); hopt = zeros(Nsym,1); hdir = hopt;
for s = 1:Nsym
    hk = H(:,1,s); gk = G(1,:,s).'; dk = D(1,1,s); cc = gk .* hk;
    hopt(s) = dk + sum(cc .* exp(1i*(angle(dk) - angle(cc))));
    hdir(s) = dk;
end
C.opt   = mean(log2(1 + Pt*abs(hopt).^2/Pn));
C.noRIS = mean(log2(1 + Pt*abs(hdir).^2/Pn));
end
