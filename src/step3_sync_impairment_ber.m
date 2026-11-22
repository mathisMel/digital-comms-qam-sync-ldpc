%% Digital Comms — Step 3: BER Degradation from CFO (ISI), Phase Offset & Timing Shift
% Course: ELEC-H401 (Digital Communications)
%
% Evaluates 4-QAM BER degradation vs. Eb/N0 across three synchronization errors:
%   1. CFO residual ISI ([0, 2, 10, 20, 50] ppm after removing phase rotation)
%   2. Uncompensated constant carrier phase offset ([0, pi/8, pi/4] rad)
%   3. Fractional symbol-timing sampling offset (t0 = [0, 0.06, 0.15, 0.30]*T)

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
Nsymb = 1e4;
RollOff = 0.2;
Rs = 5e6;
ppm = [0, 2, 10, 20, 50];   % CFO [ppm]
fc = 600e6;                 % Nominal carrier frequency [Hz]
M = 100;                    % High oversampling factor for fine timing shifts
Shift = [0, 0.02, 0.05, 0.1] * 3; % Fractional symbol timing offsets t0/T
phi0_list = [0, pi/8, pi/4];      % Constant carrier phase offsets [rad]

Fs = Rs*M;
N = 25*M;
df = (1/N)*Fs;
T = 1/Rs;
Ts = 1/Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;
SNR = 0:2:16;               % Eb/N0 sweep [dB]
CFO = ppm * 1e-6 * fc;

BERCFO = zeros(numel(ppm), numel(SNR));
BERPhase = zeros(numel(phi0_list), numel(SNR));
BERtShift = zeros(numel(Shift), numel(SNR));

%% Filter & Transmit Signal
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
Nb = Nbps * Nsymb;
inputBit = randi([0 1], [Nb, 1]);
symbTx = mapping(inputBit, Nbps, modulation);
Eb = compute_eb(symbTx, Nbps, T);
symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);

%% 1. Impact of CFO (ISI Component Only)
for j = 1:length(SNR)
    symbRxCorrupted = add_awgn(Eb, SNR(j), Rs, symbTxUpFiltered);
    t_cfo1 = ((0:length(symbRxCorrupted)-1)*Ts)';
    for p = 1:length(ppm)
        symbRxDesyncNoisy = symbRxCorrupted .* exp(1j*2*pi*CFO(p)*t_cfo1);
        symbRxCorruptedFiltered = conv(symbRxDesyncNoisy, gtMatched);
        symbRxCorruptedCropped = symbRxCorruptedFiltered(N:end-(N-1));
        t_cfo2 = ((0:length(symbRxCorruptedCropped)-1)*Ts)';
        symbRxCorruptedCroppedCFO = symbRxCorruptedCropped .* exp(-1j*2*pi*CFO(p)*t_cfo2);
        symbRxCorruptedFilteredDownCFO = downsample(symbRxCorruptedCroppedCFO, M);
        bitRxCorruptedCFO = demapping(symbRxCorruptedFilteredDownCFO, Nbps, modulation);
        BERCFO(p, j) = sum(inputBit ~= bitRxCorruptedCFO) / Nb;
    end
end

%% 2. Impact of Uncompensated Constant Phase Shift
for p_idx = 1:length(phi0_list)
    for j = 1:length(SNR)
        symbTxCorrupted = add_awgn(Eb, SNR(j), Rs, symbTxUpFiltered);
        symbRxPhase = symbTxCorrupted .* exp(1j*phi0_list(p_idx));
        symbRxFilt = conv(symbRxPhase, gtMatched);
        symbRxCrop = symbRxFilt(N:end-(N-1));
        symbRxDown = downsample(symbRxCrop, M);
        bitRxPh = demapping(symbRxDown, Nbps, modulation);
        BERPhase(p_idx, j) = sum(inputBit ~= bitRxPh) / Nb;
    end
end

%% 3. Impact of Sampling Timing Offset
for s = 1:length(Shift)
    for j = 1:length(SNR)
        symbRxCorrupted = add_awgn(Eb, SNR(j), Rs, symbTxUpFiltered);
        symbRxCorruptedFiltered = conv(symbRxCorrupted, gtMatched);
        symbRxCorruptedCropped = symbRxCorruptedFiltered(N:end-(N-1));
        shift_samples = round(Shift(s)*M);
        symbRxShifted = circshift(symbRxCorruptedCropped, shift_samples);
        symbRxCorruptedFilteredDown = downsample(symbRxShifted, M);
        bitRxCorrupted = demapping(symbRxCorruptedFilteredDown, Nbps, modulation);
        BERtShift(s, j) = sum(inputBit ~= bitRxCorrupted) / Nb;
    end
end

%% Plots
figure;
plot(real(symbTx), imag(symbTx), '*', 'DisplayName', 'Transmitted'); hold on;
plot(real(symbRxCorruptedFilteredDown), imag(symbRxCorruptedFilteredDown), 'o', 'DisplayName', 'Received (t_0 = 0.30 T, 16 dB)');
hold off; axis equal; grid on; xlabel('I'); ylabel('Q');
title('4-QAM Constellation Under 0.30 T Sampling Offset'); legend('Location', 'best');

plot_ber_curves(SNR, BERCFO, {'no CFO', 'CFO = 2 ppm', 'CFO = 10 ppm', 'CFO = 20 ppm', 'CFO = 50 ppm'}, ...
    'BER', '4-QAM BER vs. Residual CFO ISI', true);
ylim([1e-4 1]);

plot_ber_curves(SNR, BERPhase, {'\phi_0 = 0', '\phi_0 = \pi/8', '\phi_0 = \pi/4'}, ...
    'BER', '4-QAM BER vs. Constant Carrier Phase Shift', true);
ylim([1e-4 1]);

plot_ber_curves(SNR, BERtShift, {'no shift', 't_0 = 0.06 T', 't_0 = 0.15 T', 't_0 = 0.30 T'}, ...
    'BER', '4-QAM BER vs. Sampling Time Offset', true);
ylim([1e-4 1]);

fprintf('Step 3 BER @ 10 dB | BERCFO(0,50ppm): %.3e, %.3e | BERPhase(0,pi/8,pi/4): %.3e, %.3e, %.3e\n', ...
    BERCFO(1, SNR==10), BERCFO(5, SNR==10), BERPhase(:, SNR==10));

%% Export Figures (450 DPI PNG + Vector PDF)
figRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', 'figures');
pngDir = fullfile(figRoot, 'png');
pdfDir = fullfile(figRoot, 'pdf');
if ~exist(pngDir, 'dir'), mkdir(pngDir); end
if ~exist(pdfDir, 'dir'), mkdir(pdfDir); end
figs = findall(0, 'Type', 'figure');
if ~isempty(figs)
    [~, ord] = sort([figs.Number]);
    figs = figs(ord);
    for fi = 1:numel(figs)
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('step3_imp_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('step3_imp_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
