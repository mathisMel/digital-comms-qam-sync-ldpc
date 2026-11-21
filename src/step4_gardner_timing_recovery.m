%% Digital Comms — Step 4: Gardner Timing Recovery Convergence & Loop Statistics
% Course: ELEC-H401 (Digital Communications)
%
% Evaluates non-data-aided Gardner TED timing recovery under a fractional
% timing offset Shift = 0.40 T across:
%   1. Loop gain K = [0.01, 0.02, 0.03]
%   2. Carrier frequency offset ppm = [0, 10, 50]
% Produces symbol-by-symbol trajectories and 50-symbol grouped (mean +/- std)
% timing-error convergence plots.

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
ppm = [0, 10, 50];
fc = 600e6;
M = 100;
Shift = 0.4;                % Timing offset [symbols]
K = [0.01, 0.02, 0.03];     % Gardner loop gains

Fs = Rs*M;
N = 25*M;
df = (1/N)*Fs;
T = 1/Rs;
Ts = 1/Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;
CFO = ppm * 1e-6 * fc;

%% Filter & Baseband Signal
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
Nb = Nbps * Nsymb;
inputBit = randi([0 1], [Nb, 1]);
symbTx = mapping(inputBit, Nbps, modulation);
symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);
symbRx = symbTxUpFiltered;

%% 1. Sweep Loop Gain K (at CFO = 10 ppm)
epsilons = zeros(length(symbTx), length(K));
for i = 1:length(K)
    n = (1:length(symbRx)).';
    symbRxDesync = symbRx .* exp(1j*2*pi*CFO(2)*n*Ts);
    symbRxFiltered = conv(symbRxDesync, gtMatched);
    symbRxCropped = symbRxFiltered(N:end-(N-1));
    shift_samples = round(Shift * M);
    symbRxShifted = circshift(symbRxCropped, shift_samples);
    nBis = (1:length(symbRxShifted)).';
    symbRxCroppedCFO = symbRxShifted .* exp(-1j*2*pi*CFO(2)*nBis*Ts);
    [symbRxGardner, epsilons(:, i)] = gardner_ted(symbRxCroppedCFO, K(i), M);
end

%% 2. Sweep CFO ppm (at K = 0.01)
epsilonsCFO = zeros(length(symbTx), length(ppm));
for ii = 1:length(ppm)
    n = (1:length(symbRx)).';
    symbRxDesync = symbRx .* exp(1j*2*pi*CFO(ii)*n*Ts);
    symbRxFiltered = conv(symbRxDesync, gtMatched);
    symbRxCropped = symbRxFiltered(N:end-(N-1));
    shift_samples = round(Shift * M);
    symbRxShifted = circshift(symbRxCropped, shift_samples);
    nBis = (1:length(symbRxCropped)).';
    symbRxCroppedCFO = symbRxShifted .* exp(-1j*2*pi*CFO(ii)*nBis*Ts);
    [~, epsilonsCFO(:, ii)] = gardner_ted(symbRxCroppedCFO, 0.01, M);
end

%% Group Statistics (50 Symbols per Window)
groupSize = 50;
nbGroupes = size(epsilons, 1) / groupSize;
errorsOwn = T * (Shift - epsilons);
errorsCFO = T * (Shift - epsilonsCFO);

meanVec = zeros(nbGroupes, length(K));
stdvVec = zeros(nbGroupes, length(K));
meanVecCFO = zeros(nbGroupes, length(ppm));
stdvVecCFO = zeros(nbGroupes, length(ppm));
for gi = 1:nbGroupes
    sr = (gi-1)*groupSize + 1;
    er = min(gi*groupSize, size(errorsOwn, 1));
    meanVec(gi, :) = mean(errorsOwn(sr:er, :), 1);
    stdvVec(gi, :) = std(errorsOwn(sr:er, :), 1);
    meanVecCFO(gi, :) = mean(errorsCFO(sr:er, :), 1);
    stdvVecCFO(gi, :) = std(errorsCFO(sr:er, :), 1);
end

%% Plots
figure;
plot(real(symbTx), imag(symbTx), '*'); hold on;
plot(real(symbRxGardner), imag(symbRxGardner), 'o'); hold off;
axis equal; grid on; xlabel('I'); ylabel('Q');
title('4-QAM Constellation After Gardner Timing Recovery');
legend('Transmitted', 'Gardner Output (K = 0.03)', 'Location', 'best');

figure;
plot(epsilons, 'LineWidth', 1.1); hold on;
yline(Shift, 'k--', 'LineWidth', 1.2); hold off;
grid on; xlabel('Symbol Index'); ylabel('\epsilon [T]');
title('Gardner Timing Estimate vs. Loop Gain K');
legend('K = 0.01', 'K = 0.02', 'K = 0.03', 'True shift (0.4 T)', 'Location', 'southeast');

figure;
colorsOwn = lines(length(K));
for ci = 1:length(K)
    plot(meanVec(:, ci), 'o-', 'Color', colorsOwn(ci, :)); hold on;
    plot(meanVec(:, ci) + stdvVec(:, ci), '*--', 'Color', colorsOwn(ci, :));
    plot(meanVec(:, ci) - stdvVec(:, ci), '*--', 'Color', colorsOwn(ci, :));
end
hold off; grid on;
xlabel('Group Index (50 symbols/group)'); ylabel('Timing Error (mean \pm std) [s]');
title('Gardner Timing Error Convergence vs. Loop Gain K');
legend('K=0.01 mean', '+std', '-std', 'K=0.02 mean', '+std', '-std', 'K=0.03 mean', '+std', '-std', 'Location', 'northeast');

figure;
colorsCFO = lines(length(ppm));
for ci = 1:length(ppm)
    plot(meanVecCFO(:, ci), 'o-', 'Color', colorsCFO(ci, :)); hold on;
    plot(meanVecCFO(:, ci) + stdvVecCFO(:, ci), '*--', 'Color', colorsCFO(ci, :));
    plot(meanVecCFO(:, ci) - stdvVecCFO(:, ci), '*--', 'Color', colorsCFO(ci, :));
end
hold off; grid on;
xlabel('Group Index (50 symbols/group)'); ylabel('Timing Error (mean \pm std) [s]');
title('Gardner Timing Error Convergence vs. Residual CFO');
legend('0 ppm mean', '+std', '-std', '10 ppm mean', '+std', '-std', '50 ppm mean', '+std', '-std', 'Location', 'northeast');

fprintf('Step 4 Final Gardner epsilon (K=[0.01,0.02,0.03]): %s | mean err: %s s\n', ...
    sprintf('%.4f ', epsilons(end,:)), sprintf('%.2e ', meanVec(end,:)));
