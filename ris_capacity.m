function R = ris_capacity(matFile, Pt_dBm, noise_dBm)
% RIS_CAPACITY  Capacity with optimised RIS, random-phase RIS and no RIS
% from one SimRIS channel file (saved with the GUI's "Save as" button).
%
%   R = ris_capacity('InH_S1_f28_N256.mat', 30, -100)
%
% SimRIS output sizes:
%   H : N  x Nt x Nsym   (Tx -> RIS)
%   G : Nr x N  x Nsym   (RIS -> Rx)
%   D : Nr x Nt x Nsym   (direct Tx -> Rx)
% Use Nt = Nr = 1 (SISO) in the GUI for this script.

if nargin < 2, Pt_dBm = 30;     end
if nargin < 3, noise_dBm = -100; end   % same noise level as the SimRIS paper

S = load(matFile);
H = getVar(S, 'H');  G = getVar(S, 'G');  D = getVar(S, 'D');

[N, Nt, Nsym] = size(H);
Nr = size(G, 1);
if Nt ~= 1 || Nr ~= 1
    error('Run SimRIS with 1 Tx and 1 Rx antenna (SISO) for this script.');
end

Pt  = 10^((Pt_dBm   - 30)/10);   % W
Pn  = 10^((noise_dBm - 30)/10);  % W

h_opt = zeros(Nsym,1); h_rand = zeros(Nsym,1); h_dir = zeros(Nsym,1);
for k = 1:Nsym
    hk = H(:,1,k);             % N x 1
    gk = G(1,:,k).';           % N x 1
    dk = D(1,1,k);             % scalar
    c  = gk .* hk;             % cascaded per-element channel

    % Optimal phases: align every reflected path with the direct path
    theta     = angle(dk) - angle(c);
    h_opt(k)  = dk + sum(c .* exp(1j*theta));      % = |d| + sum|c| in magnitude
    % Random (unconfigured) RIS
    h_rand(k) = dk + sum(c .* exp(1j*2*pi*rand(N,1)));
    % No RIS
    h_dir(k)  = dk;
end

snr = @(h) Pt*abs(h).^2/Pn;
R.N            = N;
R.Nsym         = Nsym;
R.Pt_dBm       = Pt_dBm;
R.C_opt        = mean(log2(1 + snr(h_opt)));    % ergodic capacity, bps/Hz
R.C_rand       = mean(log2(1 + snr(h_rand)));
R.C_noRIS      = mean(log2(1 + snr(h_dir)));
R.SNR_opt_dB   = 10*log10(mean(snr(h_opt)));
R.SNR_rand_dB  = 10*log10(mean(snr(h_rand)));
R.SNR_noRIS_dB = 10*log10(mean(snr(h_dir)));
R.Prx_opt_dBm   = 10*log10(mean(Pt*abs(h_opt).^2))  + 30;
R.Prx_noRIS_dBm = 10*log10(mean(Pt*abs(h_dir).^2))  + 30;
R.RIS_gain_dB  = R.SNR_opt_dB - R.SNR_noRIS_dB;
end

function X = getVar(S, name)
f = fieldnames(S);
idx = find(strcmpi(f, name), 1);
if isempty(idx)
    error('Variable %s not found in file. Variables present: %s', name, strjoin(f', ', '));
end
X = S.(f{idx});
end
