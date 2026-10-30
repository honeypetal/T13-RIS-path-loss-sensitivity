function [H, G, D] = geo_ris_channel(env, scen, fGHz, N, Tx, Rx, RIS, Nsym, cal)
% GEO_RIS_CHANNEL  Independent 3GPP-based geometric RIS channel model
% (Team T13 - Methodology 2). Generates channels in exactly the same
% format as SimRIS so both simulators can be compared on identical
% scenarios.
%
%   [H, G, D] = geo_ris_channel('InH', 1, 28, 64, [0 25 2], [38 48 1], [40 50 2], 100)
%
%   env  : 'InH' (indoor office) or 'UMi' (urban street canyon)
%   scen : 1 = RIS on xz plane (side wall), 2 = RIS on yz plane (opposite wall)
%   fGHz : carrier frequency in GHz
%   N    : number of RIS elements (perfect square)
%   Tx, Rx, RIS : [x y z] coordinates in metres
%   Nsym : number of channel realizations
%   cal  : (optional, Methodology 3) struct with measured close-in model
%          parameters for the DIRECT link: n_LOS, sigma_LOS, n_NLOS,
%          sigma_NLOS, and optionally force_los (1 = LOS, 0 = NLOS, NaN = random)
%          Extra options (Methodology 2 ablation):
%            apply_ris = true -> use the same CI parameters on the RIS links
%            plos = '5gcm'    -> use the 5GCM LOS probability (as in SimRIS)
%
% Outputs (same sizes as SimRIS):
%   H : N x 1 x Nsym   (Tx -> RIS)
%   G : 1 x N x Nsym   (RIS -> Rx)
%   D : 1 x 1 x Nsym   (Tx -> Rx, direct)
%
% Model summary (differs from SimRIS's cluster-based model):
%   * Path loss, shadow fading and LOS probability: 3GPP TR 38.901
%     (InH-Office / UMi-Street Canyon), not the 5GCM close-in model.
%   * Small-scale fading: Rician (K = 7 dB InH, 9 dB UMi) when LOS,
%     Rayleigh when NLOS - no explicit scattering clusters.
%   * LOS phase of every RIS element computed from its exact 3-D
%     position (spherical wavefront), elements spaced lambda/2.
%   * RIS element pattern G_e*cos^(2q)(theta), G_e = pi, q = 0.285
%     (same physical element assumption as the SimRIS reference paper).
%   * Indoors the RIS links are LOS (RIS deployed in view of the AP and
%     within 3-5 m of the user); outdoors their LOS state is random.

if nargin < 9, cal = []; end
c   = 3e8;
lam = c / (fGHz*1e9);
k   = 2*pi/lam;
m   = round(sqrt(N));
if m^2 ~= N, error('N must be a perfect square'); end

% ---- RIS element positions (square array centred on RIS, spacing lam/2)
idx = ((0:m-1) - (m-1)/2) * lam/2;
[a, b] = meshgrid(idx, idx);
P = repmat(RIS(:).', N, 1);
if scen == 1          % xz plane
    P(:,1) = P(:,1) + a(:);  P(:,3) = P(:,3) + b(:);
else                  % yz plane
    P(:,2) = P(:,2) + a(:);  P(:,3) = P(:,3) + b(:);
end

if strcmpi(env, 'InH'), Kf = 10^(7/10); else, Kf = 10^(9/10); end

H = zeros(N, 1, Nsym);  G = zeros(1, N, Nsym);  D = zeros(1, 1, Nsym);

for s = 1:Nsym
    % ---- LOS states
    if strcmpi(env, 'InH')
        los1 = true;  los2 = true;                       % RIS links
    else
        los1 = rand < plos_any(env, Tx, RIS, cal);
        los2 = rand < plos_any(env, Rx, RIS, cal);
    end
    los3 = rand < plos_any(env, Tx, Rx, cal);            % direct link
    if ~isempty(cal) && isfield(cal, 'force_los') && ~isnan(cal.force_los)
        los3 = logical(cal.force_los);
    end

    H(:,1,s) = ris_link(env, Tx, P, fGHz, k, los1, Kf, cal);
    G(1,:,s) = ris_link(env, Rx, P, fGHz, k, los2, Kf, cal).';
    D(1,1,s) = direct_link(env, Tx, Rx, fGHz, k, los3, Kf, cal);
end
end

% =====================================================================
function v = ris_link(env, T, P, fGHz, k, los, Kf, cal)
% Channel between terminal T and every RIS element (rows of P).
dv = P - repmat(T(:).', size(P,1), 1);
d  = sqrt(sum(dv.^2, 2));
dc = mean(d);
hUT = min(T(3), mean(P(:,3)));          % lower end of the link
if ~isempty(cal) && isfield(cal, 'apply_ris') && cal.apply_ris && isfield(cal, 'n_NLOS')
    fs1m = 20*log10(4*pi*fGHz*1e9/3e8);
    if los, pl = fs1m + 10*cal.n_LOS*log10(dc);  sf = cal.sigma_LOS;
    else,   pl = fs1m + 10*cal.n_NLOS*log10(dc); sf = cal.sigma_NLOS; end
else
    [pl, sf] = pl_3gpp(env, los, dc, fGHz, hUT);
end
L = 10^(-(pl + randn*sf)/10);
% element radiation pattern (elevation angle of the link)
el = asind(abs(mean(dv(:,3))) / dc);
L  = L * pi * cosd(el)^(2*0.285);
n  = size(P,1);
diffuse = (randn(n,1) + 1i*randn(n,1)) / sqrt(2);
if los
    spec = exp(-1i*k*d) * exp(1i*2*pi*rand);
    v = sqrt(L) * (sqrt(Kf/(Kf+1))*spec + sqrt(1/(Kf+1))*diffuse);
else
    v = sqrt(L) * diffuse;
end
end

function h = direct_link(env, Tx, Rx, fGHz, k, los, Kf, cal)
d = norm(Tx - Rx);
if ~isempty(cal) && isfield(cal, 'n_NLOS')
    % measurement-calibrated close-in (CI) model, d0 = 1 m
    fs1m = 20*log10(4*pi*fGHz*1e9/3e8);
    if los, pl = fs1m + 10*cal.n_LOS*log10(d);  sf = cal.sigma_LOS;
    else,   pl = fs1m + 10*cal.n_NLOS*log10(d); sf = cal.sigma_NLOS; end
else
    [pl, sf] = pl_3gpp(env, los, d, fGHz, Rx(3));
end
L = 10^(-(pl + randn*sf)/10);
diffuse = (randn + 1i*randn) / sqrt(2);
if los
    h = sqrt(L) * (sqrt(Kf/(Kf+1))*exp(-1i*k*d)*exp(1i*2*pi*rand) + sqrt(1/(Kf+1))*diffuse);
else
    h = sqrt(L) * diffuse;
end
end

function [pl, sf] = pl_3gpp(env, los, d, fc, hUT)
% 3GPP TR 38.901 Table 7.4.1-1 path loss (dB) and shadow-fading std (dB)
if strcmpi(env, 'InH')
    pl_los = 32.4 + 17.3*log10(d) + 20*log10(fc);
    if los
        pl = pl_los; sf = 3.0;
    else
        pl = max(pl_los, 17.3 + 38.3*log10(d) + 24.9*log10(fc)); sf = 8.03;
    end
else  % UMi street canyon (d < breakpoint assumed)
    pl_los = 32.4 + 21*log10(d) + 20*log10(fc);
    if los
        pl = pl_los; sf = 4.0;
    else
        pl = max(pl_los, 35.3*log10(d) + 22.4 + 21.3*log10(fc) - 0.3*(hUT - 1.5)); sf = 7.82;
    end
end
end

function p = plos_any(env, A, B, cal)
% LOS probability: 3GPP (2-D distance) by default, 5GCM/SimRIS (3-D) if asked
if ~isempty(cal) && isfield(cal, 'plos') && strcmpi(cal.plos, '5gcm')
    d = norm(A - B);
    if strcmpi(env, 'InH')
        if d <= 1.2, p = 1; elseif d < 6.5, p = exp(-(d-1.2)/4.7); else, p = 0.32*exp(-(d-6.5)/32.6); end
    else
        p = min(20/d, 1)*(1 - exp(-d/39)) + exp(-d/39);
    end
else
    p = p_los(env, norm(A(1:2) - B(1:2)));
end
end

function p = p_los(env, d2)
% 3GPP TR 38.901 Table 7.4.2-1 LOS probability (d2 = 2-D distance)
if strcmpi(env, 'InH')       % InH - Mixed office
    if d2 <= 1.2
        p = 1;
    elseif d2 < 6.5
        p = exp(-(d2 - 1.2)/4.7);
    else
        p = 0.32*exp(-(d2 - 6.5)/32.6);
    end
else                         % UMi - Street canyon
    if d2 <= 18
        p = 1;
    else
        p = 18/d2 + exp(-d2/36)*(1 - 18/d2);
    end
end
end
