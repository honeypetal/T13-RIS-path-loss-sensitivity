%% Team 13 - Extra parameters from the project brief
% Completes the "Parameters to Measure" list for Project 15:
%   X1  Received power with and without RIS        (from the saved M1 datasets)
%   X2  LoS vs NLoS direct link: RIS gain           (SimRIS, direct link forced LoS / NLoS)
%   X3  Sensitivity to RIS element spacing          (SimRIS, spacing lambda/4, lambda/2, lambda)
%   X4  MIMO capacity (1x1, 2x2, 4x4)               (SimRIS, ULA at Tx and Rx)
%   X5  RMS delay spread, coherence bandwidth and angular spread (from SimRIS scatterer geometry)
%
% HOW TO RUN: put this file in C:\T13_RIS (the folder with SimRIS_GUI.m, poissrnd.m
% and the 20 SimRIS .mat files from Methodology 1) and press Run.
% Run time: roughly 10-30 minutes. Progress is printed in the Command Window.
% Results (CSV files + PNG figures) are written to the folder "results_X".
clear; clc; close all;
global SIMX
SIMX = struct();
rng(11);
outDir = 'results_X';
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;

make_simris_x();                       % callable copy of SimRIS with a small measurement hook

Pt_dBm = 30;  Pn_dBm = -100;
Pt = 10^((Pt_dBm-30)/10);  Pn = 10^((Pn_dBm-30)/10);
R  = 50;                               % realizations per point (use 100 if you have time)
c0 = 3e8;
cfg = { 'InH_S1_f28', 1, 1, 28, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S1_f73', 1, 1, 73, [0 25 2],  [38 48 1], [40 50 2];
        'InH_S2_f28', 1, 2, 28, [0 25 2],  [70 35 1], [70 30 2];
        'UMi_S1_f28', 2, 1, 28, [0 25 20], [80 75 1], [70 85 10] };
Nlist = [16 64 144 256 400];
nice = @(s) strrep(s, '_', ' ');

%% X1) Received power with and without RIS (saved M1 datasets, no new simulation)
fprintf('X1/5: received power from the saved M1 datasets...\n');
rows = cell(0, 6);
for c = 1:size(cfg,1)
    for N = Nlist
        S = load(sprintf('%s_N%d.mat', cfg{c,1}, N));
        [hdir, hopt] = siso_eff(S.H, S.G, S.D);
        rows(end+1,:) = {cfg{c,1}, N, Pt_dBm + 10*log10(mean(abs(hdir).^2)), ...
                         Pt_dBm + 10*log10(mean(abs(hopt).^2)), ...
                         10*log10(mean(abs(hopt).^2)/mean(abs(hdir).^2)), ...
                         mean(10*log10(Pt*abs(hopt).^2/Pn))}; %#ok<SAGROW>
    end
end
T1 = cell2table(rows, 'VariableNames', {'Config','N','Prx_noRIS_dBm','Prx_optRIS_dBm','Power_gain_dB','SNR_optRIS_dB'});
writetable(T1, fullfile(outDir, 'X1_received_power.csv'));  disp(T1(T1.N == 400, :));
figure('Color','w', 'Position', [100 100 900 420]); hold on; grid on;
mk = {'-o','-s','-^','-d'};
for c = 1:size(cfg,1)
    t = T1(strcmp(T1.Config, cfg{c,1}), :);
    plot(t.N, t.Prx_optRIS_dBm, mk{c}, 'LineWidth', 1.5, 'DisplayName', [nice(cfg{c,1}) ', opt. RIS']);
end
set(gca,'ColorOrderIndex',1);
for c = 1:size(cfg,1)
    t = T1(strcmp(T1.Config, cfg{c,1}), :);
    plot(t.N, t.Prx_noRIS_dBm, '--', 'LineWidth', 1.2, 'DisplayName', [nice(cfg{c,1}) ', no RIS']);
end
set(gca, 'XScale', 'log'); xlabel('Number of RIS elements N'); ylabel('Mean received power (dBm)');
title('Received power with and without RIS (P_t = 30 dBm)'); legend('Location','eastoutside','FontSize',7);
saveas(gcf, fullfile(outDir, 'X1_received_power.png'));

%% X2) RIS gain when the direct link is LoS vs NLoS
fprintf('X2/5: LoS vs NLoS direct link...\n');
rows = cell(0, 6);  states = [1 0];  stName = {'LoS','NLoS'};
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for s = 1:2
        SIMX = struct('force_los3', states(s));
        for N = Nlist
            hdir = zeros(R,1); hopt = hdir;
            for r = 1:R
                [h, g, hs] = SimRIS_x(env, scen, f, 1, N, 1, 1, Tx, Rx, RIS);
                [hdir(r), hopt(r)] = siso_eff(h(:,1), g(:,1).', hs(1,1));
            end
            Cno = mean(log2(1 + Pt*abs(hdir).^2/Pn));  Cop = mean(log2(1 + Pt*abs(hopt).^2/Pn));
            rows(end+1,:) = {name, stName{s}, N, Cno, Cop, Cop - Cno}; %#ok<SAGROW>
        end
        fprintf('   %-11s direct link %-4s: gain at N=400 = %.2f bps/Hz\n', name, stName{s}, rows{end,6});
    end
end
SIMX = struct();
T2 = cell2table(rows, 'VariableNames', {'Config','DirectLink','N','C_noRIS','C_optRIS','Gain'});
writetable(T2, fullfile(outDir, 'X2_los_vs_nlos.csv'));
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:size(cfg,1)
    subplot(2,2,c); hold on; grid on;
    a = T2(strcmp(T2.Config, cfg{c,1}) & strcmp(T2.DirectLink,'LoS'), :);
    b = T2(strcmp(T2.Config, cfg{c,1}) & strcmp(T2.DirectLink,'NLoS'), :);
    plot(a.N, a.Gain, '-o', 'LineWidth', 1.5); plot(b.N, b.Gain, '-s', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('RIS capacity gain (bps/Hz)');
    title(nice(cfg{c,1})); if c == 1, legend('Direct link LoS', 'Direct link NLoS', 'Location', 'northwest'); end
end
sgtitle('RIS capacity gain with a LoS vs an NLoS direct link');
saveas(gcf, fullfile(outDir, 'X2_los_vs_nlos.png'));

%% X3) Sensitivity to RIS element spacing (N = 144 keeps every case in SimRIS's far-field range)
fprintf('X3/5: element spacing...\n');
spac = [0.25 0.5 1.0];  Nsp = 144;  rows = cell(0, 6);
for c = [1 4]
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for sp = spac
        SIMX = struct('dis_factor', sp);
        hdir = zeros(R,1); hopt = hdir; hrnd = hdir;
        for r = 1:R
            [h, g, hs] = SimRIS_x(env, scen, f, 1, Nsp, 1, 1, Tx, Rx, RIS);
            [hdir(r), hopt(r), hrnd(r)] = siso_eff(h(:,1), g(:,1).', hs(1,1));
        end
        rows(end+1,:) = {name, sp, mean(log2(1+Pt*abs(hdir).^2/Pn)), mean(log2(1+Pt*abs(hrnd).^2/Pn)), ...
                         mean(log2(1+Pt*abs(hopt).^2/Pn)), mean(log2(1+Pt*abs(hopt).^2/Pn)) - mean(log2(1+Pt*abs(hdir).^2/Pn))}; %#ok<SAGROW>
        fprintf('   %-11s spacing %.2f lambda: opt %.2f bps/Hz\n', name, sp, rows{end,5});
    end
end
SIMX = struct();
T3 = cell2table(rows, 'VariableNames', {'Config','Spacing_lambda','C_noRIS','C_randomRIS','C_optRIS','Gain'});
writetable(T3, fullfile(outDir, 'X3_spacing.csv'));  disp(T3);
figure('Color','w', 'Position', [100 100 900 400]);
for k = 1:2
    subplot(1,2,k); c = [1 4]; t = T3(strcmp(T3.Config, cfg{c(k),1}), :);
    bar(categorical({'\lambda/4','\lambda/2','\lambda'}, {'\lambda/4','\lambda/2','\lambda'}), [t.C_noRIS t.C_randomRIS t.C_optRIS]);
    grid on; ylabel('Capacity (bps/Hz)'); title([nice(cfg{c(k),1}) ', N = 144']);
    if k == 1, legend('No RIS','Random RIS','Optimised RIS','Location','southoutside','Orientation','horizontal'); end
end
saveas(gcf, fullfile(outDir, 'X3_spacing.png'));

%% X4) MIMO capacity: 1x1, 2x2, 4x4 (ULA at Tx and Rx, equal power per antenna)
fprintf('X4/5: MIMO...\n');
Nm = [64 144 400];  Mlist = [1 2 4];  rows = cell(0, 7);
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    for M = Mlist
        for N = Nm
            Cno = zeros(R,1); Cop = Cno; Crn = Cno;
            for r = 1:R
                [h, g, hs] = SimRIS_x(env, scen, f, 1, N, M, M, Tx, Rx, RIS);  % h: N x Nt, g: N x Nr, hs: Nr x Nt
                Gt = g.';                                                     % Nr x N
                Cno(r) = mimo_cap(hs, Pt, Pn);
                Crn(r) = mimo_cap(hs + Gt*diag(exp(1i*2*pi*rand(N,1)))*h, Pt, Pn);
                Cop(r) = mimo_cap(hs + Gt*diag(exp(1i*ris_mimo_phase(h, Gt, hs)))*h, Pt, Pn);
            end
            rows(end+1,:) = {name, sprintf('%dx%d', M, M), N, mean(Cno), mean(Crn), mean(Cop), mean(Cop) - mean(Cno)}; %#ok<SAGROW>
        end
        fprintf('   %-11s %dx%d: N=400 opt %.2f bps/Hz (no RIS %.2f)\n', name, M, M, rows{end,6}, rows{end,4});
    end
end
T4 = cell2table(rows, 'VariableNames', {'Config','Antennas','N','C_noRIS','C_randomRIS','C_optRIS','Gain'});
writetable(T4, fullfile(outDir, 'X4_mimo.csv'));
figure('Color','w', 'Position', [100 100 1000 650]);
for c = 1:size(cfg,1)
    subplot(2,2,c); hold on; grid on;
    for M = Mlist
        t = T4(strcmp(T4.Config, cfg{c,1}) & strcmp(T4.Antennas, sprintf('%dx%d', M, M)), :);
        p = plot(t.N, t.C_optRIS, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('%dx%d opt. RIS', M, M));
        plot(t.N, t.C_noRIS, '--', 'Color', p.Color, 'LineWidth', 1.2, 'DisplayName', sprintf('%dx%d no RIS', M, M));
    end
    set(gca, 'XScale', 'log'); xlabel('N'); ylabel('Capacity (bps/Hz)'); title(nice(cfg{c,1}));
    if c == 1, legend('Location', 'northwest', 'FontSize', 7); end
end
sgtitle('MIMO capacity with and without RIS');
saveas(gcf, fullfile(outDir, 'X4_mimo.png'));

%% X5) RMS delay spread, coherence bandwidth and angular spread (N = 400)
fprintf('X5/5: delay and angular spread...\n');
rows = cell(0, 11);  N5 = 400;
for c = 1:size(cfg,1)
    [name, env, scen, f, Tx, Rx, RIS] = cfg{c,:};
    v = nan(R, 6);
    for r = 1:R
        SIMX = struct('capture', 1);
        [h, g, hs] = SimRIS_x(env, scen, f, 1, N5, 1, 1, Tx, Rx, RIS);
        P = SIMX.paths;                          % [length (m), power (linear), scatterer x y z]; NaN xyz = LoS
        if isempty(P), continue; end
        [dsD, asdD, asaD] = spreads(P, Tx, Rx, c0);
        Pris = sum(abs(g(:,1) .* h(:,1)))^2;     % optimised-RIS path power (coherent sum)
        Lris = norm(RIS - Tx) + norm(Rx - RIS);
        Pv = [P; Lris, Pris, RIS];               % RIS path enters as one extra tap via the RIS position
        [dsR, asdR, asaR] = spreads(Pv, Tx, Rx, c0);
        v(r,:) = [dsD, asdD, asaD, dsR, asdR, asaR];
    end
    m = median(v, 1, 'omitnan');
    rows(end+1,:) = {name, m(1)*1e9, 1/(5*m(1))/1e6, m(2), m(3), m(4)*1e9, 1/(5*m(4))/1e6, m(5), m(6), ...
                     sum(~isnan(v(:,1))), N5}; %#ok<SAGROW>
    fprintf('   %-11s DS %.1f ns -> %.1f ns with RIS\n', name, m(1)*1e9, m(4)*1e9);
end
SIMX = struct();
T5 = cell2table(rows, 'VariableNames', {'Config','DS_noRIS_ns','Bc_noRIS_MHz','ASD_noRIS_deg','ASA_noRIS_deg', ...
    'DS_optRIS_ns','Bc_optRIS_MHz','ASD_optRIS_deg','ASA_optRIS_deg','Realizations','N'});
writetable(T5, fullfile(outDir, 'X5_delay_angular_spread.csv'));  disp(T5);
figure('Color','w', 'Position', [100 100 1000 380]);
lab = categorical(cellfun(nice, cfg(:,1), 'UniformOutput', false), cellfun(nice, cfg(:,1), 'UniformOutput', false));
subplot(1,3,1); bar(lab, [T5.DS_noRIS_ns T5.DS_optRIS_ns]); grid on; ylabel('RMS delay spread (ns)'); title('Delay spread'); legend('No RIS','Opt. RIS','Location','best');
subplot(1,3,2); bar(lab, [T5.Bc_noRIS_MHz T5.Bc_optRIS_MHz]); grid on; ylabel('B_c \approx 1/(5\tau_{rms}) (MHz)'); title('Coherence bandwidth');
subplot(1,3,3); bar(lab, [T5.ASA_noRIS_deg T5.ASA_optRIS_deg]); grid on; ylabel('Azimuth spread of arrival (deg)'); title('Angular spread at Rx');
saveas(gcf, fullfile(outDir, 'X5_delay_angular_spread.png'));

fprintf('\nAll done in %.1f min. Results are in "%s".\n', toc(tStart)/60, outDir);

%% ===================== local functions =====================
function [hdir, hopt, hrnd] = siso_eff(H, G, D)
% SISO effective channel with no RIS, optimised RIS and random-phase RIS.
% H: N x 1 x K, G: 1 x N x K, D: 1 x 1 x K (or single realization)
K = size(H, 3); hdir = zeros(K,1); hopt = hdir; hrnd = hdir;
for k = 1:K
    hk = H(:,1,k); gk = G(1,:,k).'; dk = D(1,1,k); cc = gk .* hk;
    hdir(k) = dk;
    hopt(k) = dk + sum(cc .* exp(1i*(angle(dk) - angle(cc))));
    hrnd(k) = dk + sum(cc .* exp(1i*2*pi*rand(size(cc))));
end
end

function C = mimo_cap(Hm, Pt, Pn)
% Capacity with equal power over the Nt transmit antennas
[Nr, Nt] = size(Hm);
C = real(log2(det(eye(Nr) + (Pt/(Nt*Pn)) * (Hm*Hm'))));
end

function th = ris_mimo_phase(H, Gt, D)
% Alternating optimisation: align every RIS path with the dominant eigen-mode
% of the current effective channel (10 iterations).
N = size(H, 1); th = zeros(N, 1);
for it = 1:10
    Heff = D + Gt*diag(exp(1i*th))*H;
    [U, ~, V] = svd(Heff); u = U(:,1); v = V(:,1);
    a = u' * D * v;
    b = (u' * Gt).' .* (H * v);
    th = angle(a) - angle(b);
end
end

function [ds, asd, asa] = spreads(P, Tx, Rx, c0)
% RMS delay spread (s) and circular azimuth spreads (deg, 3GPP TR 38.901 Annex A)
L = P(:,1); w = P(:,2) / sum(P(:,2)); tau = L / c0;
ds = sqrt(sum(w .* tau.^2) - sum(w .* tau)^2);
xyz = P(:,3:5); los = any(isnan(xyz), 2);
dep = xyz - Tx;  arr = xyz - Rx;                 % Tx -> scatterer, Rx -> scatterer
dep(los,:) = repmat(Rx - Tx, sum(los), 1);       % LoS: Tx -> Rx
arr(los,:) = repmat(Tx - Rx, sum(los), 1);       % LoS: Rx -> Tx
phiD = atan2(dep(:,2), dep(:,1));  phiA = atan2(arr(:,2), arr(:,1));
asd = rad2deg(sqrt(-2*log(abs(sum(w .* exp(1i*phiD))))));
asa = rad2deg(sqrt(-2*log(abs(sum(w .* exp(1i*phiA))))));
end

function make_simris_x()
% Copies SimRIS's channel generator (SimRIS_v18 in SimRIS_GUI.m) into SimRIS_x.m with a
% small hook: element spacing, forced direct-link LoS/NLoS and capture of the direct-link
% paths (length, power, scatterer position). The physics of SimRIS is not changed.
txt = fileread('SimRIS_GUI.m');
i1 = strfind(txt, 'function [h,g,h_SISO]=SimRIS_v18(');
i2 = strfind(txt, 'function edit14_Callback');
if numel(i1) ~= 1 || numel(i2) ~= 1, error('Could not find SimRIS_v18 inside SimRIS_GUI.m'); end
core = txt(i1:i2-1);
nl = newline;
core = rep(core, 'function [h,g,h_SISO]=SimRIS_v18(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz)', ...
    ['function [h,g,h_SISO]=SimRIS_x(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx_xyz,RIS_xyz)' nl ...
     'global SIMX' nl 'if isempty(SIMX), SIMX = struct(); end' nl 'SIMX_P = zeros(0,5);' nl], 1);
core = rep(core, 'dis=lambda/2;', ...
    'dis=lambda/2; if isfield(SIMX,''dis_factor''), dis=SIMX.dis_factor*lambda; end', 1);
los = 'I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);';
core = rep(core, los, [los ' if isfield(SIMX,''force_los3''), I_LOS_3=SIMX.force_los3; end'], 2);
a = 'h_SISO_NLOS = h_SISO_NLOS + beta(counter2)*exp(1i*eta)*sqrt(Lcs_SISO)*transpose(array_Rx_cs_SISO(counter2,:))*array_Tx_cs_SISO(counter2,:);';
core = rep(core, a, [a ' SIMX_P(end+1,:) = [d_cs_tilde(counter2), abs(beta(counter2))^2*Lcs_SISO/M_new, Coordinates2(counter2,:)];'], 1);
b = 'h_SISO_NLOS = h_SISO_NLOS + ((randn+1i*randn)/sqrt(2))*sqrt(Lcs_SISO_3)*transpose(array_Rx_cs2(counter2,:))*array_Tx_cs2(counter2,:);';
core = rep(core, b, ['bbx = ((randn+1i*randn)/sqrt(2)); h_SISO_NLOS = h_SISO_NLOS + bbx*sqrt(Lcs_SISO_3)*transpose(array_Rx_cs2(counter2,:))*array_Tx_cs2(counter2,:);' ...
    ' SIMX_P(end+1,:) = [d_cs_3(counter2), abs(bbx)^2*Lcs_SISO_3/M_new_3, Coordinates2_3(counter2,:)];'], 1);
c1 = 'h_SISO=h_SISO_NLOS + h_SISO_LOS; % include LOS Component if I_LOS_3==1';
core = rep(core, c1, ['if I_LOS_3==1, SIMX_P(end+1,:) = [d_T_R, L_SISO_LOS, NaN, NaN, NaN]; end' nl c1 nl 'SIMX.paths = SIMX_P;'], 1);
c2 = 'h_SISO = h_SISO_NLOS + h_SISO_LOS;    % h_SISO_LOS=0  if there is no LOS';
core = rep(core, c2, ['if I_LOS_3==1, SIMX_P(end+1,:) = [d_T_R, L_SISO_LOS, NaN, NaN, NaN]; end' nl c2 nl 'SIMX.paths = SIMX_P;'], 1);
fid = fopen('SimRIS_x.m', 'w'); fwrite(fid, core); fclose(fid);
rehash;
fprintf('SimRIS_x.m created from SimRIS_GUI.m\n');
end

function s = rep(s, old, new, nExpected)
k = numel(strfind(s, old));
if k ~= nExpected
    error('SimRIS_GUI.m version not recognised: "%s" found %d times (expected %d)', old(1:min(40,end)), k, nExpected);
end
s = strrep(s, old, new);
end
