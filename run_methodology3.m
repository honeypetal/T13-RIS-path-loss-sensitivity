%% Team 13 - Methodology 3: Measurement-informed calibration
% Calibrates the DIRECT (Tx-Rx) link of both simulators - SimRIS and the
% geometric model - against published 28/73 GHz NYU WIRELESS measurements
% (Sun, MacCartney & Rappaport), reports the residual error before and
% after calibration, and re-runs all 20 RIS scenarios with the calibrated
% models.
%
% Put this file and geo_ris_channel.m (new version) in C:\T13_RIS, the
% folder that already contains SimRIS_GUI.m, poissrnd.m, the 20 SimRIS
% .mat files and the results_M2 folder. Then press Run (about 5-15 min).
clear; clc; close all;
global SIMRIS_HOOK
SIMRIS_HOOK = [];
rng(3);
outDir = 'results_M3';
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;

%% 0) Make a callable copy of SimRIS's channel generator (SimRIS_core.m)
make_simris_core();

%% 1) Published measurement reference (omnidirectional CI model, d0 = 1 m)
% Sun, MacCartney, Rappaport, "Millimeter-wave distance-dependent
% large-scale propagation measurements and path loss models for outdoor
% and indoor 5G systems", arXiv:1511.07345, Table I.
%        case     env    f   state    n    sigma  dmin  dmax
meas = {'InH28', 'InH', 28, 'LOS',  1.1,  1.8,  4.1, 21.3;
        'InH28', 'InH', 28, 'NLOS', 2.7,  9.6,  3.9, 45.9;
        'InH73', 'InH', 73, 'LOS',  1.3,  2.4,  4.1, 21.3;
        'InH73', 'InH', 73, 'NLOS', 3.2, 11.3,  3.9, 41.9;
        'UMi28', 'UMi', 28, 'LOS',  2.1,  3.5, 31.0, 54.0;
        'UMi28', 'UMi', 28, 'NLOS', 3.4,  9.7, 61.0, 186.0};
R      = 100;      % realizations per distance
nDist  = 8;        % distances per measured range
RAYL   = 10*0.5772156649/log(10);   % 2.507 dB: mean dB loss of Rayleigh fading

%% 2) Path loss of both simulators BEFORE calibration
fprintf('Step 2/5: path-loss sweep before calibration...\n');
PL = struct();
for m = 1:size(meas,1)
    [cs, env, f, st, n, sg, dmin, dmax] = meas{m,:};
    d = linspace(dmin, dmax, nDist);
    key = [cs '_' st];
    PL.(key).d = d;  PL.(key).n_meas = n;  PL.(key).s_meas = sg;
    PL.(key).sim_before = sweep('SimRIS', env, f, st, d, R, [], RAYL);
    PL.(key).geo_before = sweep('Geo',    env, f, st, d, R, [], RAYL);
    fprintf('  %-10s done (%.0f s)\n', key, toc(tStart));
end

%% 3) Calibration: fit each simulator's path-loss exponent to the measurements
fprintf('Step 3/5: fitting path-loss exponents to measurements...\n');
cases = unique(meas(:,1), 'stable');
CAL = struct();
for c = 1:numel(cases)
    cs = cases{c};
    rL = strcmp(meas(:,1), cs) & strcmp(meas(:,4), 'LOS');
    rN = strcmp(meas(:,1), cs) & strcmp(meas(:,4), 'NLOS');
    env = meas{rN,2}; f = meas{rN,3};
    % start from the measured values, measured shadowing is imposed
    p = struct('n_LOS', meas{rL,5}, 'sigma_LOS', meas{rL,6}, ...
               'n_NLOS', meas{rN,5}, 'sigma_NLOS', meas{rN,6});
    CAL.(cs).Geo = p;            % geometric model: CI structure, measured values
    % SimRIS: NLOS exponent first (it also affects LOS-forced channels), then LOS
    dN = PL.([cs '_NLOS']).d;  dL = PL.([cs '_LOS']).d;
    p.n_NLOS = secant_fit(@(x) bias_of('SimRIS', env, f, 'NLOS', dN, R, setfield(p,'n_NLOS',x), RAYL, meas{rN,5}), meas{rN,5});
    p.n_LOS  = secant_fit(@(x) bias_of('SimRIS', env, f, 'LOS',  dL, R, setfield(p,'n_LOS',x),  RAYL, meas{rL,5}), meas{rL,5});
    CAL.(cs).SimRIS = p;
    fprintf('  %s: SimRIS fitted n_LOS = %.2f, n_NLOS = %.2f (measured %.1f / %.1f)  (%.0f s)\n', ...
        cs, p.n_LOS, p.n_NLOS, meas{rL,5}, meas{rN,5}, toc(tStart));
end

%% 4) Path loss AFTER calibration (validation) and residual-error table
fprintf('Step 4/5: validating calibrated models...\n');
rows = cell(0, 9);
for m = 1:size(meas,1)
    [cs, env, f, st, n, sg] = meas{m,1:6};
    key = [cs '_' st];  d = PL.(key).d;
    PL.(key).sim_after = sweep('SimRIS', env, f, st, d, R, CAL.(cs).SimRIS, RAYL);
    PL.(key).geo_after = sweep('Geo',    env, f, st, d, R, CAL.(cs).Geo,    RAYL);
    for stage = {'before', 'after'}
        for mdl = {'sim', 'geo'}
            S = PL.(key).([mdl{1} '_' stage{1}]);
            [nf, bias, rmse, sfit] = ci_stats(S, d, f, n);
            name = 'SimRIS'; if strcmp(mdl{1}, 'geo'), name = 'Geometric'; end
            rows(end+1,:) = {key, name, stage{1}, n, nf, sg, sfit, bias, rmse}; %#ok<SAGROW>
        end
    end
end
Tpl = cell2table(rows, 'VariableNames', {'Case','Model','Stage','n_measured', ...
    'n_fitted','sigma_measured_dB','sigma_fitted_dB','Bias_dB','RMSE_dB'});
writetable(Tpl, fullfile(outDir, 'M3_pathloss_residuals.csv'));
fprintf('\n===== Direct-link path loss vs measurements =====\n'); disp(Tpl);

crow = cell(0, 6);
for c = 1:numel(cases)
    for mdl = {'SimRIS', 'Geo'}
        p = CAL.(cases{c}).(mdl{1});
        crow(end+1,:) = {cases{c}, mdl{1}, p.n_LOS, p.sigma_LOS, p.n_NLOS, p.sigma_NLOS}; %#ok<SAGROW>
    end
end
Tcal = cell2table(crow, 'VariableNames', {'Case','Model','n_LOS','sigma_LOS_dB','n_NLOS','sigma_NLOS_dB'});
writetable(Tcal, fullfile(outDir, 'M3_calibrated_parameters.csv'));
fprintf('\n===== Calibrated direct-link parameters =====\n'); disp(Tcal);

%% 5) Re-run the 20 RIS scenarios with the calibrated models
fprintf('Step 5/5: re-running the 20 RIS scenarios...\n');
Pt = 10^((30-30)/10);  Pn = 10^((-100-30)/10);
cfg = { 'InH_S1_f28', 'InH', 1, 28, [0 25 2],  [38 48 1], [40 50 2],  'InH28';
        'InH_S1_f73', 'InH', 1, 73, [0 25 2],  [38 48 1], [40 50 2],  'InH73';
        'InH_S2_f28', 'InH', 2, 28, [0 25 2],  [70 35 1], [70 30 2],  'InH28';
        'UMi_S1_f28', 'UMi', 1, 28, [0 25 20], [80 75 1], [70 85 10], 'UMi28' };
Nlist = [16 64 144 256 400];  Nsym = 100;
crows = cell(0, 14);
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS, cs] = cfg{c,:};
    envId = 1 + strcmp(env, 'UMi');
    for N = Nlist
        S0 = load(sprintf('%s_N%d.mat', name, N));                                   % SimRIS, uncalibrated (M1)
        G0 = load(fullfile('results_M2', sprintf('GEO_%s_N%d.mat', name, N)));      % Geometric, uncalibrated (M2)
        SIMRIS_HOOK = CAL.(cs).SimRIS;  SIMRIS_HOOK.force_los = NaN;
        H = zeros(N,1,Nsym); G = zeros(1,N,Nsym); D = zeros(1,1,Nsym);
        for r = 1:Nsym
            [h, g, hs] = SimRIS_core(envId, scen, f, 1, N, 1, 1, Tx, Rx, RIS);
            H(:,1,r) = h(:,1);  G(1,:,r) = g(:,1).';  D(1,1,r) = hs(1,1);
        end
        SIMRIS_HOOK = [];
        cal = CAL.(cs).Geo;  cal.force_los = NaN;
        [Hg, Gg, Dg] = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, Nsym, cal);
        a = evaluate(S0.H, S0.G, S0.D, Pt, Pn);   b = evaluate(H, G, D, Pt, Pn);
        e = evaluate(G0.Hg, G0.Gg, G0.Dg, Pt, Pn); g2 = evaluate(Hg, Gg, Dg, Pt, Pn);
        crows(end+1,:) = {name, N, a.noRIS, b.noRIS, e.noRIS, g2.noRIS, ...
                          a.opt, b.opt, e.opt, g2.opt, ...
                          a.opt-a.noRIS, b.opt-b.noRIS, e.opt-e.noRIS, g2.opt-g2.noRIS}; %#ok<SAGROW>
        fprintf('  %-11s N=%3d | no RIS: SimRIS %.2f->%.2f  Geo %.2f->%.2f | opt: SimRIS %.2f->%.2f  Geo %.2f->%.2f\n', ...
            name, N, a.noRIS, b.noRIS, e.noRIS, g2.noRIS, a.opt, b.opt, e.opt, g2.opt);
    end
end
Tc = cell2table(crows, 'VariableNames', {'Config','N', ...
    'C_noRIS_SimRIS','C_noRIS_SimRIS_cal','C_noRIS_Geo','C_noRIS_Geo_cal', ...
    'C_opt_SimRIS','C_opt_SimRIS_cal','C_opt_Geo','C_opt_Geo_cal', ...
    'Gain_SimRIS','Gain_SimRIS_cal','Gain_Geo','Gain_Geo_cal'});
writetable(Tc, fullfile(outDir, 'M3_capacity.csv'));
disp(Tc);

% cross-simulator agreement before vs after calibration
arows = cell(0, 5);
for grp = {'All', 'All_indoor'}
    sel = true(height(Tc),1);
    if strcmp(grp{1}, 'All_indoor'), sel = startsWith(Tc.Config, 'InH'); end
    for q = {'C_noRIS', 'C_opt', 'Gain'}
        rb = rmse2(Tc.([q{1} '_SimRIS'])(sel),     Tc.([q{1} '_Geo'])(sel));
        ra = rmse2(Tc.([q{1} '_SimRIS_cal'])(sel), Tc.([q{1} '_Geo_cal'])(sel));
        arows(end+1,:) = {grp{1}, q{1}, rb, ra, 100*(rb-ra)/rb}; %#ok<SAGROW>
    end
end
Ta = cell2table(arows, 'VariableNames', {'Group','Quantity','RMSE_before','RMSE_after','Reduction_pct'});
writetable(Ta, fullfile(outDir, 'M3_agreement.csv'));
fprintf('\n===== SimRIS vs geometric model: disagreement before / after calibration =====\n'); disp(Ta);

%% 6) Figures
plot_pathloss(PL, meas, f_of(meas), 'before', outDir);
plot_pathloss(PL, meas, f_of(meas), 'after',  outDir);

% fitted path-loss exponents
figure('Color','w'); keys = strcat(meas(:,1), '_', meas(:,4));
nb = zeros(numel(keys), 5);
for k = 1:numel(keys)
    sel = strcmp(Tpl.Case, keys{k});
    nb(k,:) = [meas{k,5}, ...
        Tpl.n_fitted(sel & strcmp(Tpl.Model,'SimRIS')    & strcmp(Tpl.Stage,'before')), ...
        Tpl.n_fitted(sel & strcmp(Tpl.Model,'Geometric') & strcmp(Tpl.Stage,'before')), ...
        Tpl.n_fitted(sel & strcmp(Tpl.Model,'SimRIS')    & strcmp(Tpl.Stage,'after')), ...
        Tpl.n_fitted(sel & strcmp(Tpl.Model,'Geometric') & strcmp(Tpl.Stage,'after'))];
end
bar(categorical(strrep(keys,'_',' '), strrep(keys,'_',' ')), nb); grid on;
ylabel('Effective path-loss exponent n');
legend('Measured (NYU)', 'SimRIS before', 'Geometric before', 'SimRIS after', 'Geometric after', 'Location', 'northwest');
title('Direct-link path-loss exponent: simulators vs measurements');
saveas(gcf, fullfile(outDir, 'M3_PLE_bars.png'));

% capacity per configuration
for c = 1:size(cfg,1)
    t = Tc(strcmp(Tc.Config, cfg{c,1}), :);
    figure('Color','w', 'Position', [100 100 1000 420]);
    for sp = 1:2
        subplot(1,2,sp); hold on; grid on;
        if sp == 1, q = 'C_opt'; ttl = 'Optimised RIS'; else, q = 'C_noRIS'; ttl = 'No RIS'; end
        plot(t.N, t.([q '_SimRIS']),     '-o',  'Color', [0 0.45 0.74],    'LineWidth', 1.5);
        plot(t.N, t.([q '_Geo']),        '-s',  'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
        plot(t.N, t.([q '_SimRIS_cal']), '--o', 'Color', [0 0.45 0.74],    'LineWidth', 1.5, 'MarkerFaceColor', [0 0.45 0.74]);
        plot(t.N, t.([q '_Geo_cal']),    '--s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'MarkerFaceColor', [0.85 0.33 0.10]);
        set(gca, 'XScale', 'log'); xlabel('Number of RIS elements N'); ylabel('Ergodic capacity (bps/Hz)');
        title([ttl ': ' strrep(cfg{c,1}, '_', ' ')]);
        if sp == 1, legend('SimRIS', 'Geometric', 'SimRIS calibrated', 'Geometric calibrated', 'Location', 'northwest'); end
    end
    saveas(gcf, fullfile(outDir, ['M3_Cap_' cfg{c,1} '.png']));
end

% disagreement before/after
figure('Color','w');
sel = strcmp(Ta.Group, 'All');
bar(categorical({'No-RIS capacity','Optimised RIS capacity','RIS gain'}, {'No-RIS capacity','Optimised RIS capacity','RIS gain'}), ...
    [Ta.RMSE_before(sel) Ta.RMSE_after(sel)]); grid on;
ylabel('SimRIS vs geometric RMSE (bps/Hz)'); legend('Before calibration', 'After calibration', 'Location', 'northeast');
title('Cross-simulator disagreement before and after measurement calibration');
saveas(gcf, fullfile(outDir, 'M3_agreement_bars.png'));

% RIS gain at N = 400
figure('Color','w');
t = Tc(Tc.N == 400, :);
bar(categorical(strrep(t.Config,'_',' ')), [t.Gain_SimRIS t.Gain_Geo t.Gain_SimRIS_cal t.Gain_Geo_cal]); grid on;
ylabel('Capacity gain over no-RIS (bps/Hz)');
legend('SimRIS', 'Geometric', 'SimRIS calibrated', 'Geometric calibrated', 'Location', 'northeast');
title('RIS capacity gain at N = 400 before and after calibration');
saveas(gcf, fullfile(outDir, 'M3_gain_N400.png'));

fprintf('\nDone in %.1f min. Tables and figures are in the "%s" folder.\n', toc(tStart)/60, outDir);

%% ===================== local functions =====================
function make_simris_core()
% Copies SimRIS's own channel generator (SimRIS_v18 inside SimRIS_GUI.m)
% into SimRIS_core.m and adds a hook so that the direct-link path-loss
% parameters and the direct-link LOS state can be set from outside.
txt = fileread('SimRIS_GUI.m');
i1 = strfind(txt, 'function [h,g,h_SISO]=SimRIS_v18(');
i2 = strfind(txt, 'function edit14_Callback');
if numel(i1) ~= 1 || numel(i2) ~= 1, error('Could not find SimRIS_v18 inside SimRIS_GUI.m'); end
core = txt(i1:i2-1);
core = strrep(core, 'function [h,g,h_SISO]=SimRIS_v18(', 'function [h,g,h_SISO]=SimRIS_core(');
hook = [ ...
 'global SIMRIS_HOOK' newline ...
 'if ~isempty(SIMRIS_HOOK) && isfield(SIMRIS_HOOK, ''n_NLOS'')  % T13 Methodology 3 calibration hook' newline ...
 '    if exist(''shadow'', ''var''), shadow = shadow * (SIMRIS_HOOK.sigma_NLOS / sigma_NLOS); end' newline ...
 '    n_LOS = SIMRIS_HOOK.n_LOS;   sigma_LOS = SIMRIS_HOOK.sigma_LOS;' newline ...
 '    n_NLOS = SIMRIS_HOOK.n_NLOS; sigma_NLOS = SIMRIS_HOOK.sigma_NLOS;' newline ...
 '    b_LOS = 0; b_NLOS = 0;' newline ...
 'end' newline];
if numel(strfind(core, '%% STEP 8')) ~= 1, error('SimRIS_GUI.m version not recognised (STEP 8)'); end
core = strrep(core, '%% STEP 8', ['%% STEP 8' newline hook]);
losLine = 'I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);';
force = [losLine ' if ~isempty(SIMRIS_HOOK) && isfield(SIMRIS_HOOK,''force_los'') && ~isnan(SIMRIS_HOOK.force_los), I_LOS_3 = SIMRIS_HOOK.force_los; end'];
core = strrep(core, losLine, force);
fid = fopen('SimRIS_core.m', 'w'); fwrite(fid, core); fclose(fid);
rehash;
fprintf('SimRIS_core.m created from SimRIS_GUI.m\n');
end

function PLs = sweep(model, env, f, st, d, R, cal, RAYL)
% Direct-link path loss samples (dB), R x numel(d), small-scale fading
% removed on average (Rayleigh correction for NLOS)
global SIMRIS_HOOK
if strcmp(env, 'InH'), Tx = [0 25 2];  RIS = [40 50 2];  hR = 1;
else,                  Tx = [0 25 20]; RIS = [70 85 10]; hR = 1; end
force = double(strcmp(st, 'LOS'));
PLs = zeros(R, numel(d));
for j = 1:numel(d)
    Rx = [sqrt(d(j)^2 - (Tx(3)-hR)^2), Tx(2), hR];
    if strcmp(model, 'SimRIS')
        if isempty(cal), SIMRIS_HOOK = struct(); else, SIMRIS_HOOK = cal; end
        SIMRIS_HOOK.force_los = force;
        envId = 1 + strcmp(env, 'UMi');
        for r = 1:R
            [~, ~, hs] = SimRIS_core(envId, 1, f, 1, 4, 1, 1, Tx, Rx, RIS);
            PLs(r,j) = -20*log10(abs(hs(1,1)));
        end
        SIMRIS_HOOK = [];
    else
        c = cal; if isempty(c), c = struct(); end
        c.force_los = force;
        [~, ~, D] = geo_ris_channel(env, 1, f, 1, Tx, Rx, RIS, R, c);
        PLs(:,j) = -20*log10(abs(squeeze(D)));
    end
end
if strcmp(st, 'NLOS'), PLs = PLs - RAYL; end
end

function b = bias_of(model, env, f, st, d, R, cal, RAYL, nMeas)
S = sweep(model, env, f, st, d, R, cal, RAYL);
[~, b] = ci_stats(S, d, f, nMeas);
end

function [nf, bias, rmse_d, sfit] = ci_stats(S, d, f, nMeas)
% CI fit (d0 = 1 m) of the per-distance mean path loss, plus residuals
% against the measured CI line
fs1m = 20*log10(4*pi*f*1e9/3e8);
X = 10*log10(d(:));  m = mean(S, 1).';
nf = sum((m - fs1m).*X) / sum(X.^2);
ref = fs1m + nMeas*X;
bias = mean(m - ref);  rmse_d = sqrt(mean((m - ref).^2));
res = S - repmat((fs1m + nf*X).', size(S,1), 1);
sfit = std(res(:));
end

function x = secant_fit(fun, x0)
% Solves fun(x) = 0 (zero mean bias) with a few secant steps
x1 = x0 + 0.3;  f0 = fun(x0);  f1 = fun(x1);
for it = 1:4
    if abs(f1 - f0) < 1e-6, break; end
    x2 = x1 - f1*(x1 - x0)/(f1 - f0);
    x2 = min(max(x2, 0.5), 6);
    x0 = x1; f0 = f1; x1 = x2; f1 = fun(x1);
    if abs(f1) < 0.25, break; end   % within 0.25 dB
end
x = x1;
end

function f = f_of(meas)
f = cell2mat(meas(:,3));
end

function plot_pathloss(PL, meas, fr, stage, outDir)
cases = unique(meas(:,1), 'stable');
for c = 1:numel(cases)
    figure('Color','w', 'Position', [100 100 1000 420]);
    sts = {'LOS', 'NLOS'};
    for s = 1:2
        k = find(strcmp(meas(:,1), cases{c}) & strcmp(meas(:,4), sts{s}));
        key = [cases{c} '_' sts{s}];  P = PL.(key);
        fs1m = 20*log10(4*pi*fr(k)*1e9/3e8);
        dd = linspace(min(P.d), max(P.d), 50);
        ref = fs1m + 10*P.n_meas*log10(dd);
        subplot(1,2,s); hold on; grid on;
        fill([dd fliplr(dd)], [ref+P.s_meas fliplr(ref-P.s_meas)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        plot(dd, ref, 'k-', 'LineWidth', 2);
        S1 = P.(['sim_' stage]);  S2 = P.(['geo_' stage]);
        errorbar(P.d, mean(S1,1), std(S1,0,1), 'o', 'Color', [0 0.45 0.74], 'LineWidth', 1.3, 'MarkerFaceColor', [0 0.45 0.74]);
        errorbar(P.d*1.02, mean(S2,1), std(S2,0,1), 's', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3, 'MarkerFaceColor', [0.85 0.33 0.10]);
        set(gca, 'XScale', 'log'); xlabel('Tx-Rx distance (m)'); ylabel('Path loss (dB)');
        title(sprintf('%s %s (%s calibration)', regexprep(cases{c}, '(\d+)', ' $1 GHz'), sts{s}, stage));
        legend(sprintf('Measured \\pm\\sigma (n = %.1f, \\sigma = %.1f dB)', P.n_meas, P.s_meas), 'Measured CI fit', ...
               'SimRIS (mean \pm std)', 'Geometric (mean \pm std)', 'Location', 'northwest', 'FontSize', 7);
    end
    saveas(gcf, fullfile(outDir, sprintf('M3_PL_%s_%s.png', stage, cases{c})));
end
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

function r = rmse2(a, b)
r = sqrt(mean((a(:) - b(:)).^2));
end
