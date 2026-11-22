%% Digital Comms — Step 1: Square-QAM AWGN Simulation & BER Benchmarks
% Course: ELEC-H401 (Digital Communications)
%
% Simulates baseband Gray-coded 4/16/64/256-QAM over an AWGN channel with
% Root-Raised-Cosine (RRC) pulse shaping and matched filtering (M = 4):
%   1. Generates RC/RRC impulse & frequency responses (half_nyquist_filter)
%   2. Displays 4/16/64/256-QAM constellation diagrams
%   3. Evaluates hard-decision BER across Eb/N0 = -5..25 dB

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

%% Simulation Parameters
modulation = 'qam';
Nbps = [2, 4, 6, 8];        % 4-QAM, 16-QAM, 64-QAM, 256-QAM
Nsymb = 1e4;                % Symbols per SNR point
rolloff = 0.2;              % RRC roll-off factor
Rs = 5e6;                   % Symbol rate [symbol/s]
M = 4;                      % Oversampling factor
N = 25*M;                   % Filter tap count

Fs = Rs*M;                  % Sampling frequency [Hz]
T = 1/Rs;                   % Symbol period [s]
Ts = 1/Fs;                  % Sample period [s]
df = (1/N)*Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;
SNR = -5:2:25;              % Eb/N0 sweep [dB]

BER = zeros(numel(Nbps), numel(SNR));
symbTxStorage = zeros(Nsymb, numel(Nbps));

%% Transmit / Receive RRC Filters
[gt, gtMatched] = half_nyquist_filter(T, rolloff, f, t, true);

%% Monte Carlo AWGN BER Sweep
for i = 1:length(Nbps)
    Nb = Nbps(i) * Nsymb;
    inputBit = randi([0 1], [Nb, 1]);
    symbTx = mapping(inputBit, Nbps(i), modulation);
    symbTxStorage(:, i) = symbTx;
    Eb = compute_eb(symbTx, Nbps(i), T);
    symbTxUp = upsample(symbTx, M);
    symbTxUpFiltered = conv(symbTxUp, gt);

    for j = 1:length(SNR)
        symbTxCorrupted = add_awgn(Eb, SNR(j), Rs, symbTxUpFiltered);
        symbRxCorruptedFiltered = conv(symbTxCorrupted, gtMatched);
        symbRxCorruptedCropped = symbRxCorruptedFiltered(N:end-(N-1));
        symbRxCorruptedFilteredDown = downsample(symbRxCorruptedCropped, M);
        bitRxCorrupted = demapping(symbRxCorruptedFilteredDown, Nbps(i), modulation);
        BER(i, j) = sum(inputBit ~= bitRxCorrupted) / Nb;
    end
end

%% Plots
figure;
for i = 1:size(symbTxStorage, 2)
    subplot(1, 4, i);
    plot(real(symbTxStorage(:, i)), imag(symbTxStorage(:, i)), 'o');
    xlabel('In-Phase (I)');
    if i == 1
        ylabel('Quadrature (Q)');
    end
    axis equal; grid on;
    title(['QAM ', num2str(2^(2*i))]);
end
sgtitle('Gray-Coded Square-QAM Constellations');

plot_ber_curves(SNR, BER, {'4-QAM', '16-QAM', '64-QAM', '256-QAM'}, ...
    'BER', 'Square-QAM BER over AWGN Channel', true);
ylim([1e-5 1]);

fprintf('Step 1 BER at Eb/N0 = 15 dB [4, 16, 64, 256-QAM]: %s\n', ...
    sprintf('%.3e ', BER(:, SNR == 15)));

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
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('step1_qam_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('step1_qam_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
