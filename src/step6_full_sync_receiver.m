%% Digital Comms — Step 6: Full Receiver Chain (Gardner + Pilots + Phase Interpolation)
% Course: ELEC-H401 (Digital Communications)
%
% Simulates complete 16-QAM burst reception with fractional sampling offset
% (Shift = 0.17 T) and CFO in [0, 2, 10, 50] ppm:
%   1. Periodic pilot insertion across frame
%   2. Gardner timing recovery (K = 0.01)
%   3. Frame ToA & coarse CFO acquisition via first pilot
%   4. Initial phase offset removal + inter-pilot phase unwrapping & linear interpolation

clc; clear; close all;
set(groot, 'DefaultAxesFontName', 'Calisto MT');  % Set default font for axes to Calisto MT
set(groot, 'DefaultTextFontName', 'Calisto MT');  % Set default font for text to Calisto MT
set(groot, 'DefaultFigureColor', 'w');  % Set the background color of all figures to white
try s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
rng(0);
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

%% Parameters
modulation = 'qam';
Nbps = 4;                   % 16-QAM
Nsymb = 1e4;
RollOff = 0.2;
Rs = 5e6;
ppm = [0, 2, 10, 50];
fc = 600e6;
M = 100;
Shift = 0.17;               % Sampling time shift [symbols]

Fs = Rs*M;
N = 25*M;
df = (1/N)*Fs;
T = 1/Rs;
Ts = 1/Fs;
f = -Fs/2:df:((Fs/2)-df);
t = (-N/2:(N/2)-1)*Ts;
CFO = ppm * 1e-6 * fc;

pilot_pos = 930;            % After Gardner convergence window
pilot_interval = Nsymb/100;
pilot_size = pilot_interval/5;
Number_pilots = floor((Nsymb - (floor(pilot_pos/100)*100)) / pilot_interval);
cross_correlation_avg_window = 8;

%% Filter & Pilot Frame Construction
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
Nb = Nbps * Nsymb;
inputBit = randi([0 1], [Nb, 1]);
symbTx = mapping(inputBit, Nbps, modulation);

for i = 1:Number_pilots
    startIndex = pilot_pos + (i-1)*pilot_interval;
    if i == 1
        pilotFirst = symbTx(startIndex:startIndex + pilot_size - 1);
    elseif i == 2
        pilotSecond = symbTx(startIndex:startIndex + pilot_size - 1);
    else
        symbTx(startIndex:startIndex + pilot_size - 1) = pilotSecond;
    end
end

test = pilotSecond .* conj(pilotFirst);
phase_difference_P1_P2_avg = angle(mean(test ./ abs(test)));
inputBitReal = demapping(symbTx(pilot_pos:end), Nbps, modulation);

symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);
symbRx = symbTxUpFiltered;

pilot_ToA = zeros(Number_pilots, 1);
BitErrors = zeros(numel(ppm), 1);

%% Synchronization Loop Over CFO Values
for ii = 1:length(ppm)
    t_cfo1 = ((0:length(symbRx)-1)*Ts)';
    symbRxDesync = symbRx .* exp(1j*2*pi*CFO(ii)*t_cfo1);
    symbRxFiltered = conv(symbRxDesync, gtMatched);
    symbRxCropped = symbRxFiltered(N:end-(N-1));
    shift_samples = round(Shift * M);
    symbRxShifted = circshift(symbRxCropped, shift_samples);

    [symbRxGardnerCFO, epsilon] = gardner_ted(symbRxShifted, 0.01, M);

    pilotPhasesDiff = zeros(pilot_size, Number_pilots);
    ReceivedPilot = zeros(pilot_size, Number_pilots);

    [toa, est_CFO] = frame_cfo_acquisition(symbRxGardnerCFO, pilotFirst, cross_correlation_avg_window, T);
    symbRxNoISICFO = symbRxGardnerCFO(toa:end);
    if ii > 1
        t_cfo2 = ((0:length(symbRxNoISICFO)-1)*Ts)';
        symbRxNoISICFO = symbRxNoISICFO .* exp((-1j*2*pi)*est_CFO*t_cfo2);
    end

    toa_new = 1;
    ReceivedPilot(:, 1) = symbRxNoISICFO(toa_new:toa_new + pilot_size - 1);
    initial_phase_offset_CFO = mean(angle(ReceivedPilot(:, 1) .* conj(pilotFirst)), 1);
    symbRxNoISICFO = symbRxNoISICFO .* exp(-1j*initial_phase_offset_CFO);

    ReceivedPilot(:, 1) = symbRxNoISICFO(toa_new:toa_new + pilot_size - 1);
    pilotPhasesDiff(:, 1) = 0;
    pilot_ToA(1, :) = toa_new;

    for i = 2:Number_pilots
        idxStart = toa_new + (i-1)*pilot_interval;
        ReceivedPilot(:, i) = symbRxNoISICFO(idxStart:idxStart + pilot_size - 1);
        pilotPhasesDiff(:, i) = angle(ReceivedPilot(:, i) .* conj(ReceivedPilot(:, i-1)));
        pilot_ToA(i, :) = idxStart;
    end

    pilotPhasesDiff = angle(mean(exp(1j*pilotPhasesDiff), 1));
    pilotPhasesDiff(1, 2) = pilotPhasesDiff(1, 2) - phase_difference_P1_P2_avg;
    pilotPhasesDiff = pilotPhasesDiff';

    % Unwrap negative inter-pilot phase increments into [0, 2*pi) for monotonic CFO drift
    if CFO(ii) ~= 0
        for i = 1:length(pilotPhasesDiff)
            if pilotPhasesDiff(i, 1) < 0
                pilotPhasesDiff(i, 1) = pi + (pi - abs(pilotPhasesDiff(i, 1)));
            end
        end
    end
    pilotPhasesDiffUnwrapped = cumsum(pilotPhasesDiff);

    dataTimes = (1:length(symbRxNoISICFO)).';
    interpolatedPhases = interp1(pilot_ToA, pilotPhasesDiffUnwrapped, dataTimes, 'linear', 'extrap');
    symbRxNoISICFO = symbRxNoISICFO .* exp(-1j*interpolatedPhases);

    bitRxCorruptedCFO = demapping(symbRxNoISICFO, Nbps, modulation);
    BitErrors(ii, 1) = sum(inputBitReal ~= bitRxCorruptedCFO);

    if ii == 1 || ii == 3
        figure;
        plot(epsilon); hold on;
        yline(Shift, 'r--', 'LineWidth', 1.5); hold off;
        ylabel('Timing offset \epsilon [T]'); xlabel('Symbol Index'); grid on;
        title(sprintf('Gardner Timing Lock (%d ppm CFO)', ppm(ii)));

        figure;
        plot(real(symbTx), imag(symbTx), '*', 'DisplayName', 'Transmitted'); hold on;
        plot(real(symbRxNoISICFO), imag(symbRxNoISICFO), 'o', 'DisplayName', 'Equalized'); hold off;
        axis equal; grid on; xlabel('I'); ylabel('Q');
        title(sprintf('Synchronized 16-QAM Constellation (%d ppm CFO)', ppm(ii)));
        legend('Location', 'best');
    end
end

fprintf('Step 6 Summary | Final toa = %d | est_CFO = %.2f Hz | Bit errors [0,2,10,50 ppm]: %s\n', ...
    toa, est_CFO, sprintf('%d ', BitErrors));
