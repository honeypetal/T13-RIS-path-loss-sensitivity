%% Team 13 - Methodology 2 extension, Step B: statistical confidence
% Part 1 - Convergence: how many channel realizations are needed before
%          the capacity estimate stops changing? (geometric model, N = 400)
% Part 2 - Bootstrap 95% confidence intervals for every SimRIS and
%          geometric result, and a test of which SimRIS-vs-geometric
%          differences are larger than Monte Carlo noise.
% Run inside C:\T13_RIS (needs the SimRIS .mat files and results_M2).
clear; clc; close all;
outDir = 'results_M2X';  if ~exist(outDir, 'dir'), mkdir(outDir); end
Pt = 1; Pn = 10^((-100-30)/10);
rng(7);

cfg = { 'InH_S1_f28', 'InH', 1, 28, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S1_f73', 'InH', 1, 73, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S2_f28', 'InH', 2, 28, [0 25 2],  [70 35 1], [70 30 2];
        'UMi_S1_f28', 'UMi', 1, 28, [0 25 20], [80 75 1], [70 85 10] };

%% Part 1: convergence of the Monte Carlo estimate
NsymList = [25 50 100 200 500 1000];
trials   = 20;          % independent repetitions for each Nsym
N = 400;
conv = cell(0, 7);
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for ns = NsymList
        co = zeros(trials,1); cn = co;
        for t = 1:trials
            [H, G, D] = geo_ris_channel(env, scen, f, N, Tx, Rx, RIS, ns);
            [r, ~] = evaluate(H, G, D, Pt, Pn);
            co(t) = r.opt; cn(t) = r.noRIS;
        end
        conv(end+1,:) = {name, ns, mean(co), std(co), mean(cn), std(cn), std(co - cn)}; %#ok<SAGROW>
    end
    fprintf('Convergence %-11s done\n', name);
end
Tconv = cell2table(conv, 'VariableNames', {'Config','Nsym','C_opt_mean','C_opt_std', ...
    'C_noRIS_mean','C_noRIS_std','Gain_std'});
writetable(Tconv, fullfile(outDir, 'B_convergence.csv'));
disp(Tconv);

figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:size(cfg,1)
    t = Tconv(strcmp(Tconv.Config, cfg{c,1}), :);
    subplot(2,2,c); hold on; grid on;
    shade(t.Nsym, t.C_opt_mean,   1.96*t.C_opt_std,   [0 0.45 0.74]);
    shade(t.Nsym, t.C_noRIS_mean, 1.96*t.C_noRIS_std, [0.85 0.33 0.10]);
    xline(100, 'k:', 'Nsym = 100', 'LabelVerticalAlignment', 'bottom');
    set(gca, 'XScale', 'log'); xlabel('Number of channel realizations'); ylabel('Capacity (bps/Hz)');
    title([strrep(cfg{c,1}, '_', ' ') ', N = 400']);
    if c == 1, legend('Opt RIS 95% band', 'Opt RIS mean', 'No RIS 95% band', 'No RIS mean', 'Location', 'best', 'FontSize', 7); end
end
saveas(gcf, fullfile(outDir, 'B_convergence.png'));

figure('Color','w'); hold on; grid on;
for c = 1:size(cfg,1)
    t = Tconv(strcmp(Tconv.Config, cfg{c,1}), :);
    loglog(t.Nsym, t.C_opt_std, '-o', 'LineWidth', 1.5, 'DisplayName', strrep(cfg{c,1}, '_', ' '));
end
x = NsymList; y = Tconv.C_opt_std(1) * sqrt(NsymList(1) ./ x);
loglog(x, y, 'k--', 'DisplayName', '1/\surd{N_{sym}} reference');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of channel realizations'); ylabel('Std. dev. of C_{opt} estimate (bps/Hz)');
title('Monte Carlo error decreases as 1/\surd{N_{sym}}'); legend('Location', 'northeast');
saveas(gcf, fullfile(outDir, 'B_convergence_std.png'));

%% Part 2: bootstrap 95% confidence intervals and significance test
Bn = 1000;   Nlist = [16 64 144 256 400];
rows = cell(0, 14);
for c = 1:size(cfg,1)
    name = cfg{c,1};
    for N = Nlist
        S  = load(sprintf('%s_N%d.mat', name, N));
        Gm = load(fullfile('results_M2', sprintf('GEO_%s_N%d.mat', name, N)));
        [~, xs] = evaluate(S.H, S.G, S.D, Pt, Pn);
        [~, xg] = evaluate(Gm.Hg, Gm.Gg, Gm.Dg, Pt, Pn);
        [ms_o, ls_o, us_o] = boot(xs.c_opt, Bn);   [mg_o, lg_o, ug_o] = boot(xg.c_opt, Bn);
        [ms_n, ls_n, us_n] = boot(xs.c_no, Bn);    [mg_n, lg_n, ug_n] = boot(xg.c_no, Bn);
        % bootstrap CI of the difference (Geo - SimRIS)
        [~, ld_o, ud_o] = bootdiff(xg.c_opt, xs.c_opt, Bn);
        [~, ld_n, ud_n] = bootdiff(xg.c_no,  xs.c_no,  Bn);
        rows(end+1,:) = {name, N, ms_o, ls_o, us_o, mg_o, lg_o, ug_o, ...
                         ~(ld_o <= 0 && ud_o >= 0), ms_n, mg_n, ld_n, ud_n, ~(ld_n <= 0 && ud_n >= 0)}; %#ok<SAGROW>
    end
end
Tb = cell2table(rows, 'VariableNames', {'Config','N','C_opt_SimRIS','C_opt_SimRIS_lo','C_opt_SimRIS_hi', ...
    'C_opt_Geo','C_opt_Geo_lo','C_opt_Geo_hi','C_opt_diff_significant', ...
    'C_noRIS_SimRIS','C_noRIS_Geo','C_noRIS_diff_lo','C_noRIS_diff_hi','C_noRIS_diff_significant'});
writetable(Tb, fullfile(outDir, 'B_bootstrap.csv'));
disp(Tb);
fprintf('\nOptimised-RIS capacity: %d of %d SimRIS-vs-geometric differences are statistically significant (95%%).\n', ...
    sum(Tb.C_opt_diff_significant), height(Tb));
fprintf('No-RIS capacity:        %d of %d differences are statistically significant (95%%).\n', ...
    sum(Tb.C_noRIS_diff_significant), height(Tb));

figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:size(cfg,1)
    t = Tb(strcmp(Tb.Config, cfg{c,1}), :);
    subplot(2,2,c); hold on; grid on;
    errorbar(t.N, t.C_opt_SimRIS, t.C_opt_SimRIS - t.C_opt_SimRIS_lo, t.C_opt_SimRIS_hi - t.C_opt_SimRIS, ...
        '-o', 'LineWidth', 1.4, 'Color', [0 0.45 0.74]);
    errorbar(t.N*1.04, t.C_opt_Geo, t.C_opt_Geo - t.C_opt_Geo_lo, t.C_opt_Geo_hi - t.C_opt_Geo, ...
        '--s', 'LineWidth', 1.4, 'Color', [0.30 0.75 0.93]);
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('Optimised RIS capacity (bps/Hz)');
    title([strrep(cfg{c,1}, '_', ' ') ' (95% bootstrap CI)']);
    if c == 1, legend('SimRIS', 'Geometric', 'Location', 'northwest'); end
end
saveas(gcf, fullfile(outDir, 'B_bootstrap_CI.png'));
fprintf('\nStep B done. Results in "%s".\n', outDir);

%% ---- local functions ----
function [C, x] = evaluate(H, G, D, Pt, Pn)
Nsym = size(H, 3); hopt = zeros(Nsym,1); hdir = hopt;
for s = 1:Nsym
    hk = H(:,1,s); gk = G(1,:,s).'; dk = D(1,1,s); cc = gk .* hk;
    hopt(s) = dk + sum(cc .* exp(1i*(angle(dk) - angle(cc))));
    hdir(s) = dk;
end
x.c_opt = log2(1 + Pt*abs(hopt).^2/Pn);
x.c_no  = log2(1 + Pt*abs(hdir).^2/Pn);
C.opt = mean(x.c_opt);  C.noRIS = mean(x.c_no);
end

function [m, lo, hi] = boot(x, B)
n = numel(x); bm = zeros(B,1);
for b = 1:B, bm(b) = mean(x(randi(n, n, 1))); end
bm = sort(bm); m = mean(x);
lo = bm(max(1, round(0.025*B)));  hi = bm(round(0.975*B));
end

function [m, lo, hi] = bootdiff(x, y, B)
nx = numel(x); ny = numel(y); bm = zeros(B,1);
for b = 1:B, bm(b) = mean(x(randi(nx, nx, 1))) - mean(y(randi(ny, ny, 1))); end
bm = sort(bm); m = mean(x) - mean(y);
lo = bm(max(1, round(0.025*B)));  hi = bm(round(0.975*B));
end

function shade(x, m, h, col)
x = x(:).'; m = m(:).'; h = h(:).';
fill([x fliplr(x)], [m+h fliplr(m-h)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, m, '-o', 'Color', col, 'LineWidth', 1.5);
end
