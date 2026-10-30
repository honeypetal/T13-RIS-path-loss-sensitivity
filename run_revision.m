%% Team 13 - Revision run (answers the second-level corrections)
%  A  Hook check: SimRIS with our hook (inactive) gives exactly the same channels as unmodified SimRIS
%  B  Direct-link path-loss diagnostics: per-realization path loss of SimRIS vs its own law and NYU,
%     street-canyon (UMi) path-length check (SimRIS uses a_c + |Tx - scatterer| instead of
%     a_c + |scatterer - Rx|), LoS leakage check, effective exponent with and without the correction
%  C  Main study with COMMON RANDOM NUMBERS (same direct channel for every N and every parameter set),
%     K = 1000 realizations, three direct-link parameter sets (NYU / 5GCM / 3GPP exponents) in both
%     models, ablation for all four configurations, equivalence tests (+-0.5 bps/Hz), bandwidth check
%  D  Sensitivity sweep: RIS gain vs direct-link NLOS exponent
%  E  Direct link forced LoS / NLoS (same seeds)
%
% HOW TO RUN: put this file in C:\T13_RIS (with SimRIS_GUI.m, poissrnd.m and geo_ris_channel.m from
% Methodology 3) and press Run. First try QUICK = true (about 3-5 min) to check everything works,
% then set QUICK = false for the full run (roughly 40-90 min). Results go to the folder "results_REV".
clear; clc; close all;
QUICK = true;                % <-- set to false for the full run
RUN_A = true; RUN_B = true; RUN_C = true; RUN_D = true; RUN_E = true;

global SIMX
SIMX = struct();
outDir = 'results_REV'; if ~exist(outDir, 'dir'), mkdir(outDir); end
tAll = tic;
make_cores();

if QUICK, K = 40; Ksw = 20; Klos = 20; Kd = 40; Kb = 200; else, K = 1000; Ksw = 500; Klos = 500; Kd = 1000; Kb = 1000; end
Pt = 1;                                   % 30 dBm
noise_dBm = [-100 -87 -81];  BWlab = {'5 MHz','100 MHz','400 MHz'};   % -174 dBm/Hz + 10log10(B) + NF 7 dB
Pn = 10.^((noise_dBm - 30)/10);
cfg = { 'InH_S1_f28', 'InH', 1, 28, [0 25 2],  [38 48 1], [40 50 2], 256;
        'InH_S1_f73', 'InH', 1, 73, [0 25 2],  [38 48 1], [40 50 2], 400;
        'InH_S2_f28', 'InH', 2, 28, [0 25 2],  [70 35 1], [70 30 2], 400;
        'UMi_S1_f28', 'UMi', 1, 28, [0 25 20], [80 75 1], [70 85 10], 400 };   % last column: headline N (far field)
Nlist = [16 64 144 256 400];
% direct-link parameter sets (close-in form, FSPL(1 m) anchor):  [n_LOS sigma_LOS n_NLOS sigma_NLOS]
setName = {'NYU', '5GCM', '3GPP'};
PS.InH28 = [1.1 1.8 2.7 9.6;   1.73 3.02 3.19*(1+0.06*(28-24.2)/24.2) 8.29;  1.73 3.0 3.83 8.03];
PS.InH73 = [1.3 2.4 3.2 11.3;  1.73 3.02 3.19*(1+0.06*(73-24.2)/24.2) 8.29;  1.73 3.0 3.83 8.03];
PS.UMi28 = [2.1 3.5 3.4 9.7;   1.98 3.1  3.19 8.2;                           2.1  4.0 3.53 7.82];
psKey = @(env, f) sprintf('%s%d', env, f);
seed = @(c, r) 100000*c + r;
fs1m = @(f) 20*log10(4*pi*f*1e9/3e8);

%% ===================== A) hook check =====================
if RUN_A
fprintf('A) Hook check...\n');
rows = cell(0, 5);
for c = [1 4]
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c, 1:7}; envId = 1 + strcmp(env, 'UMi');
    for r = 1:10
        SIMX = struct();
        rng(seed(c, r)); [h1, g1, d1] = SimRIS_ref(envId, scen, f, 1, 64, 1, 1, Tx, Rx, RIS);
        rng(seed(c, r)); [h2, g2, d2] = SimRIS_r(envId, scen, f, 1, 64, 1, 1, Tx, Rx, RIS);
        rows(end+1,:) = {name, r, max(abs(h1(:)-h2(:))), max(abs(g1(:)-g2(:))), max(abs(d1(:)-d2(:)))}; %#ok<SAGROW>
    end
end
TA = cell2table(rows, 'VariableNames', {'Config','Seed','MaxDiff_H','MaxDiff_G','MaxDiff_D'});
writetable(TA, fullfile(outDir, 'A_hook_check.csv'));
fprintf('   largest difference between hooked and unmodified SimRIS: %g\n', max([TA.MaxDiff_H; TA.MaxDiff_G; TA.MaxDiff_D]));
end

%% ===================== B) direct-link path-loss diagnostics =====================
if RUN_B
fprintf('B) Direct-link path-loss diagnostics...\n');
dgs = { 'InH28', 'InH', 28, [0 25 2],  3.9, 45.9, 3.19*(1+0.06*(28-24.2)/24.2), 2.7;
         'UMi28', 'UMi', 28, [0 25 20], 61,  186,  3.19,                         3.4 };
RAYL = 10*0.5772156649/log(10);         % mean dB loss of Rayleigh fading (2.507 dB)
raw = cell(0, 5); sm = cell(0, 13);
for q = 1:size(dgs, 1)
    [cs, env, f, Tx, dmin, dmax, nLaw, nNYU] = dgs{q,:}; envId = 1 + strcmp(env, 'UMi');
    dlist = logspace(log10(dmin), log10(dmax), 8);
    modes = {'native'}; if strcmp(env, 'UMi'), modes = {'native', 'corrected'}; end
    for mo = 1:numel(modes)
        PLmean = zeros(size(dlist));
        for j = 1:numel(dlist)
            d = dlist(j); Rx = [sqrt(d^2 - (Tx(3)-1)^2), Tx(2), 1];
            if strcmp(env, 'InH'), RIS = [40 50 2]; else, RIS = [70 85 10]; end
            pl = zeros(Kb,1); nLOSrows = 0; ratio = []; shortFrac = [];
            for r = 1:Kb
                SIMX = struct('force_los3', 0, 'capture', 1, 'fix_umi', strcmp(modes{mo}, 'corrected'));
                rng(seed(10+q, r));
                [~, ~, hs] = SimRIS_r(envId, 1, f, 1, 16, 1, 1, Tx, Rx, RIS);
                pl(r) = -20*log10(abs(hs(1,1)));
                P = SIMX.paths;
                if ~isempty(P)
                    nLOSrows = nLOSrows + sum(any(isnan(P(:,3:5)), 2));
                    S = P(~any(isnan(P(:,3:5)), 2), :);
                    geo = vecnorm(S(:,3:5) - Tx, 2, 2) + vecnorm(S(:,3:5) - Rx, 2, 2);
                    ratio = [ratio; S(:,1) ./ geo]; %#ok<AGROW>
                    shortFrac = [shortFrac; S(:,1) < d]; %#ok<AGROW>
                end
                raw(end+1,:) = {cs, modes{mo}, d, r, pl(r)}; %#ok<SAGROW>
            end
            PLmean(j) = mean(pl) - RAYL;
            law = fs1m(f) + 10*nLaw*log10(d);  nyu = fs1m(f) + 10*nNYU*log10(d);
            sm(end+1,:) = {cs, modes{mo}, d, law, nyu, mean(pl) - RAYL, median(pl), -10*log10(mean(10.^(-pl/10))), ...
                           std(pl), nLOSrows, mean(ratio), mean(shortFrac), NaN}; %#ok<SAGROW>
        end
        neff = sum((PLmean - fs1m(f)) .* (10*log10(dlist))) / sum((10*log10(dlist)).^2);   % CI fit through FSPL(1 m)
        for k2 = size(sm,1)-numel(dlist)+1:size(sm,1), sm{k2, 13} = neff; end
        fprintf('   %s %-9s effective NLOS exponent %.2f (law %.2f, NYU %.2f)\n', cs, modes{mo}, neff, nLaw, nNYU);
    end
end
TB = cell2table(sm, 'VariableNames', {'Case','Mode','Distance_m','Law_dB','NYU_dB','MeanPL_dB_RaylCorr','MedianPL_dB', ...
    'PowerMeanPL_dB','StdPL_dB','LoS_paths_found','MeanPathRatio_SimRIS_to_geometric','FracPathsShorterThanTxRx','EffectiveExponent'});
writetable(TB, fullfile(outDir, 'B_pathloss_summary.csv'));
TBr = cell2table(raw, 'VariableNames', {'Case','Mode','Distance_m','Realization','PL_dB'});
writetable(TBr, fullfile(outDir, 'B_pathloss_per_realization.csv'));
figure('Color','w','Position',[100 100 1100 420]);
for q = 1:2
    subplot(1,2,q); hold on; grid on;
    cs = dgs{q,1}; t = TBr(strcmp(TBr.Case, cs), :); ms = unique(t.Mode, 'stable');
    for mo = 1:numel(ms)
        u = t(strcmp(t.Mode, ms{mo}), :);
        resid = u.PL_dB - (fs1m(dgs{q,3}) + 10*dgs{q,7}*log10(u.Distance_m));
        [ys, ~] = sort(resid); plot(ys, (1:numel(ys))/numel(ys), 'LineWidth', 1.5, 'DisplayName', ['SimRIS ' ms{mo}]);
    end
    xline(0, 'k--', 'DisplayName', 'SimRIS own law');
    xlabel('Per-realization path loss minus SimRIS law (dB)'); ylabel('CDF'); title([cs ' NLOS direct link']);
    legend('Location', 'southeast');
end
saveas(gcf, fullfile(outDir, 'B_pathloss_cdf.png'));
end

%% ===================== C) main study, common random numbers =====================
if RUN_C
fprintf('C) Main study with common random numbers (K = %d)...\n', K);
rows = cell(0, 14); dcheck = cell(0, 4); bw = cell(0, 8); store = struct();
tC = tic; nCalls = 0;
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS, Nh] = cfg{c,:}; envId = 1 + strcmp(env, 'UMi'); P = PS.(psKey(env, f));
    variants = {'SimRIS', 'Geo'};
    if strcmp(env, 'UMi'), variants = {'SimRIS', 'Geo', 'SimRIS_corr'}; end
    for s = 1:3
        par = P(s,:); calG = struct('n_LOS', par(1), 'sigma_LOS', par(2), 'n_NLOS', par(3), 'sigma_NLOS', par(4), 'force_los', NaN);
        for v = 1:numel(variants)
            if strcmp(variants{v}, 'SimRIS_corr') && s ~= 2, continue; end
            Dref = [];
            for N = Nlist
                hd = zeros(K,1); ho = hd; hr = hd;
                for r = 1:K
                    switch variants{v}
                        case {'SimRIS', 'SimRIS_corr'}
                            SIMX = struct('n_LOS', par(1), 'sigma_LOS', par(2), 'n_NLOS', par(3), 'sigma_NLOS', par(4), ...
                                          'fix_umi', strcmp(variants{v}, 'SimRIS_corr'));
                            rng(seed(c, r)); [h, g, d] = SimRIS_r(envId, scen, f, 1, N, 1, 1, Tx, Rx, RIS);
                            h = h(:,1); g = g(:,1); d = d(1,1); nCalls = nCalls + 1;
                        case 'Geo'
                            rng(seed(c, r));       [~, ~, D1] = geo_ris_channel(env, scen, f, 1, Tx, Rx, RIS, 1, calG);
                            rng(seed(c, r) + 5e6); [H1, G1]   = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, 1, calG);
                            h = H1(:,1,1); g = G1(1,:,1).'; d = D1(1,1,1);
                    end
                    cc = g .* h;
                    hd(r) = d;  ho(r) = d + sum(abs(cc)) * exp(1i*angle(d));
                    rng(seed(c, r) + 9e6); hr(r) = d + sum(cc .* exp(1i*2*pi*rand(N,1)));
                end
                if isempty(Dref), Dref = hd; else, dcheck(end+1,:) = {name, variants{v}, setName{s}, max(abs(hd - Dref))}; end %#ok<SAGROW>
                Cn = log2(1 + Pt*abs(hd).^2/Pn(1)); Co = log2(1 + Pt*abs(ho).^2/Pn(1)); Cr = log2(1 + Pt*abs(hr).^2/Pn(1));
                gci = bootci_mean(Co - Cn, 1000, 0.95);
                key = sprintf('%s_%s_%s_N%d', name, variants{v}, setName{s}, N);
                store.(matlab.lang.makeValidName(key)) = struct('Cn', Cn, 'Co', Co, 'Cr', Cr);
                rows(end+1,:) = {name, variants{v}, setName{s}, N, N == Nh, mean(Cn), mean(Cr), mean(Co), mean(Co - Cn), gci(1), gci(2), ...
                                 mean(Cr - Cn), 100*mean(10*log10(Pt*abs(hd).^2/Pn(1)) < 10), 100*mean(10*log10(Pt*abs(ho).^2/Pn(1)) < 10)}; %#ok<SAGROW>
                if N == Nh
                    for b = 1:3
                        bw(end+1,:) = {name, variants{v}, setName{s}, BWlab{b}, noise_dBm(b), mean(log2(1 + Pt*abs(hd).^2/Pn(b))), ...
                                       mean(log2(1 + Pt*abs(ho).^2/Pn(b))), mean(log2(1 + Pt*abs(ho).^2/Pn(b)) - log2(1 + Pt*abs(hd).^2/Pn(b)))}; %#ok<SAGROW>
                    end
                end
            end
            fprintf('   %-11s %-11s %-4s done  (%.1f min so far)\n', name, variants{v}, setName{s}, toc(tC)/60);
        end
    end
    % ablation variants of the geometric model against SimRIS with the 5GCM set
    par = P(2,:);
    abl = {'V0_3GPP_native', [], 'V1_5GCM_pathloss', struct('n_LOS',par(1),'sigma_LOS',par(2),'n_NLOS',par(3),'sigma_NLOS',par(4),'force_los',NaN), ...
           'V2_5GCM_pathloss_and_PLOS', struct('n_LOS',par(1),'sigma_LOS',par(2),'n_NLOS',par(3),'sigma_NLOS',par(4),'force_los',NaN,'plos','5gcm')};
    for a = 1:2:numel(abl)
        if strcmp(abl{a}, 'V1_5GCM_pathloss'), continue; end    % identical to the Geo / 5GCM run above
        for N = Nlist
            hd = zeros(K,1); ho = hd;
            for r = 1:K
                rng(seed(c, r));       [~, ~, D1] = geo_ris_channel(env, scen, f, 1, Tx, Rx, RIS, 1, abl{a+1});
                rng(seed(c, r) + 5e6); [H1, G1]   = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, 1, abl{a+1});
                cc = G1(1,:,1).' .* H1(:,1,1); hd(r) = D1; ho(r) = D1 + sum(abs(cc))*exp(1i*angle(D1));
            end
            Cn = log2(1 + Pt*abs(hd).^2/Pn(1)); Co = log2(1 + Pt*abs(ho).^2/Pn(1));
            rows(end+1,:) = {name, ['Geo_' abl{a}], '5GCM', N, N == Nh, mean(Cn), NaN, mean(Co), mean(Co - Cn), NaN, NaN, NaN, NaN, NaN}; %#ok<SAGROW>
        end
    end
end
TC = cell2table(rows, 'VariableNames', {'Config','Model','ParamSet','N','Headline','C_noRIS','C_randomRIS','C_optRIS','Gain', ...
    'Gain_CI95_lo','Gain_CI95_hi','RandomRIS_minus_noRIS','Outage_noRIS_pct','Outage_optRIS_pct'});
writetable(TC, fullfile(outDir, 'C_main_summary.csv'));
writetable(cell2table(dcheck, 'VariableNames', {'Config','Model','ParamSet','MaxChange_of_D_across_N'}), fullfile(outDir, 'C_common_random_numbers_check.csv'));
writetable(cell2table(bw, 'VariableNames', {'Config','Model','ParamSet','Bandwidth','Noise_dBm','C_noRIS','C_optRIS','Gain'}), fullfile(outDir, 'C_bandwidth.csv'));
save(fullfile(outDir, 'C_per_realization.mat'), 'store', '-v7.3');

% equivalence tests (SimRIS vs Geo, same parameter set): 90% bootstrap CI of the difference inside +-0.5 bps/Hz
eq = cell(0, 9);
for c = 1:size(cfg,1)
    for s = 1:3
        for N = Nlist
            a = store.(matlab.lang.makeValidName(sprintf('%s_SimRIS_%s_N%d', cfg{c,1}, setName{s}, N)));
            b = store.(matlab.lang.makeValidName(sprintf('%s_Geo_%s_N%d', cfg{c,1}, setName{s}, N)));
            qn = {'C_noRIS', a.Cn, b.Cn; 'C_optRIS', a.Co, b.Co; 'Gain', a.Co - a.Cn, b.Co - b.Cn};
            for k2 = 1:3
                ci = bootci_diff(qn{k2,2}, qn{k2,3}, 1000, 0.90);
                eq(end+1,:) = {cfg{c,1}, setName{s}, N, qn{k2,1}, mean(qn{k2,2}) - mean(qn{k2,3}), ci(1), ci(2), ...
                               ci(1) > -0.5 && ci(2) < 0.5, N == cfg{c,8}}; %#ok<SAGROW>
            end
        end
    end
end
TE = cell2table(eq, 'VariableNames', {'Config','ParamSet','N','Quantity','Diff_SimRIS_minus_Geo','CI90_lo','CI90_hi','Equivalent_within_0p5','Headline'});
writetable(TE, fullfile(outDir, 'C_equivalence.csv'));

% ablation summary for ALL four configurations (RMSE over N against SimRIS, 5GCM set)
ab = cell(0, 5); vlist = {'Geo_V0_3GPP_native', 'Geo', 'Geo_V2_5GCM_pathloss_and_PLOS'};
vlab = {'V0 3GPP (native)', 'V1 + 5GCM path loss', 'V2 + 5GCM LOS probability'};
for c = 0:size(cfg,1)
    if c == 0, sel = true(height(TC),1); cname = 'All'; else, sel = strcmp(TC.Config, cfg{c,1}); cname = cfg{c,1}; end
    ref = TC(sel & strcmp(TC.Model, 'SimRIS') & strcmp(TC.ParamSet, '5GCM'), :);
    for v = 1:3
        t = TC(sel & strcmp(TC.Model, vlist{v}) & strcmp(TC.ParamSet, '5GCM'), :);
        ab(end+1,:) = {cname, vlab{v}, rmse(t.C_noRIS, ref.C_noRIS), rmse(t.C_optRIS, ref.C_optRIS), rmse(t.Gain, ref.Gain)}; %#ok<SAGROW>
    end
end
writetable(cell2table(ab, 'VariableNames', {'Config','Variant','RMSE_noRIS','RMSE_optRIS','RMSE_gain'}), fullfile(outDir, 'C_ablation_all_configs.csv'));

% figure: RIS gain vs N, three parameter sets, both models
figure('Color','w','Position',[100 100 1100 700]); col = lines(3);
for c = 1:size(cfg,1)
    subplot(2,2,c); hold on; grid on;
    for s = 1:3
        a = TC(strcmp(TC.Config,cfg{c,1}) & strcmp(TC.Model,'SimRIS') & strcmp(TC.ParamSet,setName{s}), :);
        b = TC(strcmp(TC.Config,cfg{c,1}) & strcmp(TC.Model,'Geo') & strcmp(TC.ParamSet,setName{s}), :);
        plot(a.N, a.Gain, '-o', 'Color', col(s,:), 'LineWidth', 1.5, 'DisplayName', ['SimRIS, ' setName{s}]);
        plot(b.N, b.Gain, '--s', 'Color', col(s,:), 'LineWidth', 1.2, 'DisplayName', ['Geometric, ' setName{s}]);
    end
    if c == 1, xline(256, ':', 'far-field limit', 'HandleVisibility', 'off'); end
    set(gca,'XScale','log'); xlabel('N'); ylabel('RIS capacity gain (bps/Hz)'); title(strrep(cfg{c,1},'_',' '));
    if c == 1, legend('Location','northwest','FontSize',7); end
end
saveas(gcf, fullfile(outDir, 'C_gain_vs_N_three_sets.png'));
fprintf('   main study done: %d SimRIS calls in %.1f min\n', nCalls, toc(tC)/60);
end

%% ===================== D) sensitivity sweep over the NLOS exponent =====================
if RUN_D
fprintf('D) Sensitivity sweep (K = %d)...\n', Ksw);
nsw = [2.4 2.7 3.0 3.2 3.4 3.6 3.83 4.0]; rows = cell(0, 8);
for c = 1:3
    [name, env, scen, f, Tx, Rx, RIS, Nh] = cfg{c,:}; envId = 1; P = PS.(psKey(env, f));
    for n = nsw
        for v = 1:2
            hd = zeros(Ksw,1); ho = hd;
            for r = 1:Ksw
                if v == 1
                    SIMX = struct('n_LOS', P(2,1), 'sigma_LOS', P(2,2), 'n_NLOS', n, 'sigma_NLOS', P(2,4));
                    rng(seed(c, r)); [h, g, d] = SimRIS_r(envId, scen, f, 1, Nh, 1, 1, Tx, Rx, RIS); h = h(:,1); g = g(:,1); d = d(1,1);
                else
                    calG = struct('n_LOS', P(2,1), 'sigma_LOS', P(2,2), 'n_NLOS', n, 'sigma_NLOS', P(2,4), 'force_los', NaN);
                    rng(seed(c, r));       [~, ~, D1] = geo_ris_channel(env, scen, f, 1, Tx, Rx, RIS, 1, calG);
                    rng(seed(c, r) + 5e6); [H1, G1]   = geo_ris_channel(env, scen, f, Nh, Tx, Rx, RIS, 1, calG);
                    h = H1(:,1,1); g = G1(1,:,1).'; d = D1(1,1,1);
                end
                hd(r) = d; ho(r) = d + sum(abs(g .* h)) * exp(1i*angle(d));
            end
            Cn = log2(1 + Pt*abs(hd).^2/Pn(1)); Co = log2(1 + Pt*abs(ho).^2/Pn(1)); ci = bootci_mean(Co - Cn, 1000, 0.95);
            mdl = {'SimRIS','Geo'};
            rows(end+1,:) = {name, Nh, n, mdl{v}, mean(Cn), mean(Co), mean(Co - Cn), (ci(2)-ci(1))/2}; %#ok<SAGROW>
        end
    end
    fprintf('   %s done\n', name);
end
TD = cell2table(rows, 'VariableNames', {'Config','N','n_NLOS','Model','C_noRIS','C_optRIS','Gain','Gain_CI95_halfwidth'});
writetable(TD, fullfile(outDir, 'D_gain_vs_exponent.csv'));
figure('Color','w','Position',[100 100 1100 380]);
for c = 1:3
    subplot(1,3,c); hold on; grid on;
    for v = 1:2
        mdl = {'SimRIS','Geo'}; t = TD(strcmp(TD.Config, cfg{c,1}) & strcmp(TD.Model, mdl{v}), :);
        errorbar(t.n_NLOS, t.Gain, t.Gain_CI95_halfwidth, '-o', 'LineWidth', 1.4, 'DisplayName', strrep(mdl{v},'Geo','Geometric'));
    end
    P = PS.(psKey(cfg{c,2}, cfg{c,4}));
    for s = 1:3, xline(P(s,3), ':', setName{s}, 'HandleVisibility', 'off'); end
    xlabel('Direct-link NLOS exponent n'); ylabel('RIS capacity gain (bps/Hz)');
    title(sprintf('%s, N = %d', strrep(cfg{c,1},'_',' '), cfg{c,8})); if c == 1, legend('Location','northwest'); end
end
saveas(gcf, fullfile(outDir, 'D_gain_vs_exponent.png'));
end

%% ===================== E) direct link forced LoS / NLoS =====================
if RUN_E
fprintf('E) Forced LoS / NLoS direct link (K = %d)...\n', Klos);
rows = cell(0, 7); st = [1 0]; stN = {'LoS','NLoS'};
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS, Nh] = cfg{c,:}; envId = 1 + strcmp(env, 'UMi'); P = PS.(psKey(env, f));
    for k2 = 1:2
        hd = zeros(Klos,1); ho = hd;
        for r = 1:Klos
            SIMX = struct('n_LOS', P(2,1), 'sigma_LOS', P(2,2), 'n_NLOS', P(2,3), 'sigma_NLOS', P(2,4), 'force_los3', st(k2));
            rng(seed(c, r)); [h, g, d] = SimRIS_r(envId, scen, f, 1, Nh, 1, 1, Tx, Rx, RIS);
            hd(r) = d(1,1); ho(r) = d(1,1) + sum(abs(g(:,1).*h(:,1))) * exp(1i*angle(d(1,1)));
        end
        Cn = log2(1 + Pt*abs(hd).^2/Pn(1)); Co = log2(1 + Pt*abs(ho).^2/Pn(1));
        rows(end+1,:) = {name, Nh, stN{k2}, mean(Cn), mean(Co), mean(Co - Cn), 100*(mean(Co)-mean(Cn))/mean(Cn)}; %#ok<SAGROW>
    end
end
writetable(cell2table(rows, 'VariableNames', {'Config','N','DirectLink','C_noRIS','C_optRIS','Gain','Gain_pct'}), fullfile(outDir, 'E_los_vs_nlos.csv'));
end
SIMX = struct();
fprintf('\nAll done in %.1f min. Results are in "%s".\n', toc(tAll)/60, outDir);

%% ===================== local functions =====================
function ci = bootci_mean(x, B, lev)
x = x(:); n = numel(x); m = zeros(B,1);
for b = 1:B, m(b) = mean(x(randi(n, n, 1))); end
m = sort(m); a = (1 - lev)/2; ci = [m(max(1, round(a*B))), m(min(B, round((1-a)*B)))];
end

function ci = bootci_diff(x, y, B, lev)
x = x(:); y = y(:); nx = numel(x); ny = numel(y); m = zeros(B,1);
for b = 1:B, m(b) = mean(x(randi(nx, nx, 1))) - mean(y(randi(ny, ny, 1))); end
m = sort(m); a = (1 - lev)/2; ci = [m(max(1, round(a*B))), m(min(B, round((1-a)*B)))];
end

function e = rmse(a, b)
e = sqrt(mean((a(:) - b(:)).^2));
end

function make_cores()
% SimRIS_ref.m : unmodified copy of SimRIS's channel generator (SimRIS_v18 in SimRIS_GUI.m)
% SimRIS_r.m   : same code plus switches that are OFF unless a field of the global SIMX is set:
%                direct-link close-in parameters, forced direct-link LoS/NLoS, capture of the
%                direct-link paths, and an optional correction of the street-canyon path length.
txt = fileread('SimRIS_GUI.m');
i1 = strfind(txt, 'function [h,g,h_SISO]=SimRIS_v18(');
i2 = strfind(txt, 'function edit14_Callback');
if numel(i1) ~= 1 || numel(i2) ~= 1, error('Could not find SimRIS_v18 inside SimRIS_GUI.m'); end
core = txt(i1:i2-1);
sig = 'function [h,g,h_SISO]=SimRIS_v18(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz)';
ref = rep(core, sig, strrep(sig, 'SimRIS_v18', 'SimRIS_ref'), 1);
fid = fopen('SimRIS_ref.m', 'w'); fwrite(fid, ref); fclose(fid);
nl = newline;
r = rep(core, sig, [strrep(sig, 'SimRIS_v18', 'SimRIS_r') nl 'global SIMX' nl 'if isempty(SIMX), SIMX = struct(); end' nl 'SIMX_P = zeros(0,5);' nl], 1);
hook = ['if isfield(SIMX, ''n_NLOS'')' nl ...
        '    if exist(''shadow'', ''var''), shadow = shadow * (SIMX.sigma_NLOS / sigma_NLOS); end' nl ...
        '    n_LOS = SIMX.n_LOS;   sigma_LOS = SIMX.sigma_LOS;' nl ...
        '    n_NLOS = SIMX.n_NLOS; sigma_NLOS = SIMX.sigma_NLOS;' nl ...
        '    b_LOS = 0; b_NLOS = 0;' nl 'end' nl];
r = rep(r, '%% STEP 8', ['%% STEP 8' nl hook], 1);
los = 'I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);';
r = rep(r, los, [los ' if isfield(SIMX,''force_los3''), I_LOS_3=SIMX.force_los3; end'], 2);
bl = 'b_cs_3(counter2)=norm(Tx_xyz-Coordinates2_3(counter2,:));';
r = rep(r, bl, ['if isfield(SIMX,''fix_umi'') && SIMX.fix_umi, b_cs_3(counter2)=norm(Rx_xyz-Coordinates2_3(counter2,:)); else, ' bl ' end'], 1);
a = 'h_SISO_NLOS = h_SISO_NLOS + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);';
r = rep(r, a, [a ' SIMX_P(end+1,:) = [d_cs_tilde(counter2), abs(beta(counter2))^2*Lcs_SISO/M_new, Coordinates2(counter2,:)];'], 1);
b = 'h_SISO_NLOS = h_SISO_NLOS + ((randn+1i*randn)/sqrt(2))*sqrt(Lcs_SISO_3)*transpose(array_Rx_cs2(counter2,:))*array_Tx_cs2(counter2,:);';
r = rep(r, b, ['bbx = ((randn+1i*randn)/sqrt(2)); h_SISO_NLOS = h_SISO_NLOS + bbx*sqrt(Lcs_SISO_3)*transpose(array_Rx_cs2(counter2,:))*array_Tx_cs2(counter2,:);' ...
    ' SIMX_P(end+1,:) = [d_cs_3(counter2), abs(bbx)^2*Lcs_SISO_3/M_new_3, Coordinates2_3(counter2,:)];'], 1);
c1 = 'h_SISO=h_SISO_NLOS + h_SISO_LOS; % include LOS Component if I_LOS_3==1';
r = rep(r, c1, ['if I_LOS_3==1, SIMX_P(end+1,:) = [d_T_R, L_SISO_LOS, NaN, NaN, NaN]; end' nl c1 nl 'SIMX.paths = SIMX_P;'], 1);
c2 = 'h_SISO = h_SISO_NLOS + h_SISO_LOS;    % h_SISO_LOS=0  if there is no LOS';
r = rep(r, c2, ['if I_LOS_3==1, SIMX_P(end+1,:) = [d_T_R, L_SISO_LOS, NaN, NaN, NaN]; end' nl c2 nl 'SIMX.paths = SIMX_P;'], 1);
fid = fopen('SimRIS_r.m', 'w'); fwrite(fid, r); fclose(fid);
rehash;
fprintf('SimRIS_ref.m and SimRIS_r.m created from SimRIS_GUI.m\n');
end

function s = rep(s, old, new, nExpected)
k = numel(strfind(s, old));
if k ~= nExpected
    error('SimRIS_GUI.m version not recognised: "%s" found %d times (expected %d)', old(1:min(40,end)), k, nExpected);
end
s = strrep(s, old, new);
end
