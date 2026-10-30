%% Team 13 - Methodology 1: RIS vs no-RIS capacity using SimRIS
% Put this script, ris_capacity.m and all your saved SimRIS .mat files
% in ONE folder, then press Run.
%
% File naming rule (important, the script reads it):
%   <Env>_<Scenario>_f<GHz>_N<elements>.mat
%   e.g.  InH_S1_f28_N64.mat   UMi_S1_f73_N1024.mat
clear; clc; close all;

Pt_dBm    = 30;      % transmit power (change if you like)
noise_dBm = -100;    % noise power, as used in the SimRIS paper
outDir    = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end

files = dir('*_N*.mat');
if isempty(files), error('No SimRIS .mat files found in this folder.'); end

%% 1) Evaluate every file -> results table
rows = cell(0, 12);
for i = 1:numel(files)
    fn  = files(i).name;
    tok = regexp(fn, '^(\w+?)_(S\d)_f(\d+)_N(\d+)\.mat$', 'tokens', 'once');
    if isempty(tok)
        warning('Skipping %s (name does not follow the rule)', fn); continue;
    end
    R = ris_capacity(fn, Pt_dBm, noise_dBm);
    rows(end+1,:) = {tok{1}, tok{2}, str2double(tok{3}), R.N, ...
        R.C_noRIS, R.C_rand, R.C_opt, R.SNR_noRIS_dB, R.SNR_opt_dB, ...
        R.Prx_noRIS_dBm, R.Prx_opt_dBm, R.RIS_gain_dB}; %#ok<SAGROW>
    fprintf('%-26s  noRIS %.2f | random %.2f | optimised %.2f bps/Hz\n', ...
        fn, R.C_noRIS, R.C_rand, R.C_opt);
end
T = cell2table(rows, 'VariableNames', {'Env','Scenario','Freq_GHz','N', ...
    'C_noRIS','C_randomRIS','C_optRIS','SNR_noRIS_dB','SNR_optRIS_dB', ...
    'Prx_noRIS_dBm','Prx_optRIS_dBm','RIS_gain_dB'});
T.C_gain    = T.C_optRIS - T.C_noRIS;                 % bps/Hz
T.C_gain_pct = 100 * T.C_gain ./ T.C_noRIS;          % percent
T = sortrows(T, {'Env','Scenario','Freq_GHz','N'});
disp(T);
writetable(T, fullfile(outDir, 'M1_results.csv'));

%% 2) Plot: capacity vs number of RIS elements (one figure per setup)
setups = unique(strcat(T.Env, '_', T.Scenario, '_f', string(T.Freq_GHz)));
for s = 1:numel(setups)
    m = strcat(T.Env, '_', T.Scenario, '_f', string(T.Freq_GHz)) == setups(s);
    t = T(m, :);
    figure('Name', setups(s));
    plot(t.N, t.C_optRIS, '-o', 'LineWidth', 1.5); hold on;
    plot(t.N, t.C_randomRIS, '-s', 'LineWidth', 1.5);
    plot(t.N, t.C_noRIS, '--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log'); grid on;
    xlabel('Number of RIS elements N'); ylabel('Ergodic capacity (bps/Hz)');
    title(strrep(setups(s), '_', ' '));
    legend('RIS (optimised phases)', 'RIS (random phases)', 'No RIS', 'Location', 'northwest');
    saveas(gcf, fullfile(outDir, "Cap_vs_N_" + setups(s) + ".png"));
end

%% 3) Plot: capacity vs transmit power for the largest-N file of each setup
PtSweep = -10:5:40;
for s = 1:numel(setups)
    m = strcat(T.Env, '_', T.Scenario, '_f', string(T.Freq_GHz)) == setups(s);
    t = T(m, :);
    fn = sprintf('%s_%s_f%d_N%d.mat', t.Env{end}, t.Scenario{end}, t.Freq_GHz(end), t.N(end));
    Copt = zeros(size(PtSweep)); Cno = Copt;
    for p = 1:numel(PtSweep)
        R = ris_capacity(fn, PtSweep(p), noise_dBm);
        Copt(p) = R.C_opt; Cno(p) = R.C_noRIS;
    end
    figure('Name', "Pt sweep " + setups(s));
    plot(PtSweep, Copt, '-o', PtSweep, Cno, '--', 'LineWidth', 1.5); grid on;
    xlabel('Transmit power P_t (dBm)'); ylabel('Ergodic capacity (bps/Hz)');
    title(sprintf('%s, N = %d', strrep(setups(s), '_', ' '), t.N(end)));
    legend('RIS (optimised)', 'No RIS', 'Location', 'northwest');
    saveas(gcf, fullfile(outDir, "Cap_vs_Pt_" + setups(s) + ".png"));
end

%% 4) Plot: capacity gain (optimised RIS - no RIS) vs N for every setup
figure('Name', 'Capacity gain comparison'); hold on; grid on;
mk = {'-o', '-s', '-^', '-d', '-v', '-p'};
for s = 1:numel(setups)
    m = strcat(T.Env, '_', T.Scenario, '_f', string(T.Freq_GHz)) == setups(s);
    t = T(m, :);
    plot(t.N, t.C_gain, mk{mod(s-1, numel(mk)) + 1}, 'LineWidth', 1.5, ...
        'DisplayName', strrep(setups(s), '_', ' ') + " GHz");
end
set(gca, 'XScale', 'log');
xlabel('Number of RIS elements N'); ylabel('Capacity gain over no-RIS (bps/Hz)');
legend('Location', 'northwest'); title('RIS capacity gain vs N for all setups');
saveas(gcf, fullfile(outDir, 'Gain_vs_N_all.png'));

%% 5) Plot: 28 GHz vs 73 GHz, same geometry (InH, Scenario 1)
m28 = strcmp(T.Env, 'InH') & strcmp(T.Scenario, 'S1') & T.Freq_GHz == 28;
m73 = strcmp(T.Env, 'InH') & strcmp(T.Scenario, 'S1') & T.Freq_GHz == 73;
if any(m28) && any(m73)
    figure('Name', 'Frequency comparison'); hold on; grid on;
    t = T(m28, :); plot(t.N, t.C_gain, '-o', 'LineWidth', 1.5, 'DisplayName', '28 GHz');
    t = T(m73, :); plot(t.N, t.C_gain, '-s', 'LineWidth', 1.5, 'DisplayName', '73 GHz');
    set(gca, 'XScale', 'log');
    xlabel('Number of RIS elements N'); ylabel('Capacity gain over no-RIS (bps/Hz)');
    legend('Location', 'northwest'); title('RIS capacity gain vs N: 28 vs 73 GHz (InH, Scenario 1)');
    saveas(gcf, fullfile(outDir, 'Gain_vs_N_freq.png'));
end

fprintf('\nDone. Table and figures are in the "%s" folder.\n', outDir);
