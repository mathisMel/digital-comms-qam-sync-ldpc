%% Digital Comms — Step 2: Phase & CFO Constellation Rotation Effects
% Course: ELEC-H401 (Digital Communications)
%
% Visualizes non-noisy 4-QAM constellation distortions under:
%   1. Constant carrier phase offset (phi_0 = pi/8 rad)
%   2. Carrier frequency offset (CFO) causing progressive phase rotation across symbols

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

%% Parameters
modulation = 'qam';
BasicBitSequence = [1; 1; 0; 1; 0; 0; 1; 0];
Nbps = 2;                   % 4-QAM
Nsymb = 1e3;
RollOff = 0.2;
Rs = 5e6;
ppm = 1;                    % CFO [ppm]
fc = 600e6;                 % Carrier frequency [Hz]
N = 100;
M = 4;
phi_0 = pi/8;               % Constant phase error [rad]

Fs = Rs*M;
df = (1/N)*Fs;
T = 1/Rs;
Ts = 1/Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;
CFO = ppm * 1e-6 * fc;

%% Pulse Shaping & Baseband Transmission
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
Nb = Nbps * Nsymb;
Nrep = Nb / length(BasicBitSequence);
inputBits = repmat(BasicBitSequence, Nrep, 1);
symbTx = mapping(inputBits, Nbps, modulation);
symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);
symbRx = symbTxUpFiltered;

%% Case A: Constant Carrier Phase Shift (phi_0 = pi/8)
symbRxPhaseFilt = conv(symbRx .* exp(1j*phi_0), gtMatched);
symbRxPhaseDown = downsample(symbRxPhaseFilt(N:end-(N-1)), M);

figure;
plot(real(symbTx), imag(symbTx), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(real(symbRxPhaseDown), imag(symbRxPhaseDown), '*', 'MarkerSize', 7);
hold off; axis equal; grid on; xlim([-1 1]); ylim([-1 1]);
xlabel('Real Part (I)'); ylabel('Imaginary Part (Q)');
title('4-QAM Constellation Under Constant Phase Offset (\phi_0 = \pi/8)');
legend('Transmitted', 'Received (shifted)', 'Location', 'best');

%% Case B: Progressive CFO Phase Drift
n = (1:length(symbRx)).';
symbRxDesync = symbRx .* exp(1j*2*pi*CFO*n*T);
symbRxCorruptedFiltered = conv(symbRxDesync, gtMatched);
symbRxCorruptedCropped = symbRxCorruptedFiltered(N:end-(N-1));
symbRxCorruptedFilteredDown = downsample(symbRxCorruptedCropped, M);

figure;
plot(real(symbTx), imag(symbTx), 'o', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
numShown = min(40, length(symbRxCorruptedFilteredDown));
scatter(real(symbRxCorruptedFilteredDown(1:numShown)), ...
        imag(symbRxCorruptedFilteredDown(1:numShown)), ...
        45, 1:numShown, 'filled');
colorbar;
hold off; axis equal; grid on; xlim([-1 1]); ylim([-1 1]);
xlabel('Real Part (I)'); ylabel('Imaginary Part (Q)');
title('Progressive Constellation Rotation Under CFO (First 40 Symbols)');
legend('Reference symbols', 'Rotating symbols (index)', 'Location', 'best');

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
        exportgraphics(figs(fi), fullfile(pngDir, sprintf('step2_rot_fig%02d.png', fi)), 'Resolution', 450, 'BackgroundColor', 'white');
        exportgraphics(figs(fi), fullfile(pdfDir, sprintf('step2_rot_fig%02d.pdf', fi)), 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
