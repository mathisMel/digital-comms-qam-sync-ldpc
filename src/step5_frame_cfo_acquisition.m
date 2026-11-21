%% Digital Comms — Step 5: Frame (ToA) & CFO Acquisition Performance vs. N and K
% Course: ELEC-H401 (Digital Communications)
%
% Evaluates differential cross-correlation frame acquisition (ToA) and CFO
% estimation standard deviation across Eb/N0 = 0..16 dB for:
%   1. Pilot symbol lengths N = [10, 20, 40] symbols (at fixed window K = 8)
%   2. Averaging windows K = [1, 8, 16] (at fixed pilot length N = 20 symbols)

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

%% Parameters
modulation = 'qam';
Nbps = 2;                   % 4-QAM
Nsymb = 1e3;
RollOff = 0.2;
Rs = 5e6;
fc = 600e6;
M = 100;
Nexp = 20;                  % Monte Carlo trials per SNR point

SNR = 0:2:16;               % Eb/N0 sweep [dB]
Fs = Rs*M;
N = 25*M;
df = (1/N)*Fs;
T = 1/Rs;
Ts = 1/Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;

pilotPos = 34;                              % True start symbol of pilot
pilotSizeSymb = [10, 20, 40];               % Pilot lengths [symbols] (N)
crossCorrelationAvgWindow = [1, 8, 16];     % Differential lag windows (K < N)
cfo_true = 10e-6 * fc;                      % 10 ppm CFO for CFO std test

%% Filter & Baseband Frame
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
Nb = Nbps * Nsymb;
inputBit = randi([0 1], [Nb, 1]);

pilotCell = cell(length(pilotSizeSymb), 1);
for i = 1:length(pilotSizeSymb)
    pilotBits = inputBit((pilotPos-1)*Nbps + (1:pilotSizeSymb(i)*Nbps));
    pilotCell{i} = mapping(pilotBits, Nbps, modulation);
end

symbTx = mapping(inputBit, Nbps, modulation);
Eb = compute_eb(symbTx, Nbps, T);
symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);

t_cfo = ((0:length(symbTxUpFiltered)-1)*Ts)';
symbCFO_base = symbTxUpFiltered .* exp(1j*2*pi*cfo_true*t_cfo);

%% 1. ToA Estimation Error (Mean |ToA_hat - ToA|)
ToAErrPilot = zeros(numel(pilotSizeSymb), numel(SNR));
ToAErrWindow = zeros(numel(crossCorrelationAvgWindow), numel(SNR));
cfoErrPilot = zeros(numel(pilotSizeSymb), numel(SNR), Nexp);
cfoErrWindow = zeros(numel(crossCorrelationAvgWindow), numel(SNR), Nexp);

for e = 1:Nexp
    for i = 1:length(pilotSizeSymb)
        for j = 1:length(SNR)
            rxN = add_awgn(Eb, SNR(j), Rs, symbCFO_base);
            rxF = conv(rxN, gtMatched);
            rxD = downsample(rxF(N:end-(N-1)), M);
            [toa, est_c] = frame_cfo_acquisition(rxD, pilotCell{i}, crossCorrelationAvgWindow(2), T);
            ToAErrPilot(i, j) = ToAErrPilot(i, j) + abs(toa - pilotPos);
            cfoErrPilot(i, j, e) = (est_c - cfo_true) / fc * 1e6;
        end
    end

    for i = 1:length(crossCorrelationAvgWindow)
        for j = 1:length(SNR)
            rxN = add_awgn(Eb, SNR(j), Rs, symbCFO_base);
            rxF = conv(rxN, gtMatched);
            rxD = downsample(rxF(N:end-(N-1)), M);
            [toa, est_c] = frame_cfo_acquisition(rxD, pilotCell{2}, crossCorrelationAvgWindow(i), T);
            ToAErrWindow(i, j) = ToAErrWindow(i, j) + abs(toa - pilotPos);
            cfoErrWindow(i, j, e) = (est_c - cfo_true) / fc * 1e6;
        end
    end
end

ToAErrPilot = ToAErrPilot ./ Nexp;
ToAErrWindow = ToAErrWindow ./ Nexp;
cfoStdPilot = std(cfoErrPilot, 0, 3);
cfoStdWindow = std(cfoErrWindow, 0, 3);

%% Plots
plot_ber_curves(SNR, ToAErrPilot, {'N = 10', 'N = 20', 'N = 40'}, ...
    'Mean |ToA Error| [symbols]', 'Frame Synchronization Error vs. Pilot Length (K = 8)', false);
plot_ber_curves(SNR, ToAErrWindow, {'K = 1', 'K = 8', 'K = 16'}, ...
    'Mean |ToA Error| [symbols]', 'Frame Synchronization Error vs. Averaging Window K (N = 20)', false);

plot_ber_curves(SNR, cfoStdPilot, {'N = 10', 'N = 20', 'N = 40'}, ...
    'CFO Estimation Std [ppm]', 'CFO Estimation Std vs. Pilot Length (K = 8)', true);
plot_ber_curves(SNR, cfoStdWindow, {'K = 1', 'K = 8', 'K = 16'}, ...
    'CFO Estimation Std [ppm]', 'CFO Estimation Std vs. Averaging Window K (N = 20)', true);

fprintf('Step 5 @ Eb/N0 = 16 dB | ToA err N=[10,20,40]: %s | CFO std [ppm]: %s\n', ...
    sprintf('%.2f ', ToAErrPilot(:, end)), sprintf('%.2f ', cfoStdPilot(:, end)));
