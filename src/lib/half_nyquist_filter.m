function [gt, gtMatched] = half_nyquist_filter(T, roll_off, f, t, doPlot)
%% Root-Raised-Cosine (Half-Nyquist) Pulse-Shaping Filter Generator
% Course: ELEC-H401 (Digital Communications)
%
% Builds a discrete Raised-Cosine spectrum H(f), normalizes so the time-domain
% RC peak h(0) = 1, takes G(f) = sqrt(H(f)), and transforms to the RRC impulse
% response g(t) as a column vector (with matched filter flipud(gt)).
%
% Inputs:
%   T        - Symbol period [s]
%   roll_off - Roll-off factor beta in [0, 1]
%   f        - Frequency grid [Hz]
%   t        - Time grid [s] (for optional plotting)
%   doPlot   - Optional boolean (default true) to display filter response figure
%
% Outputs:
%   gt        - RRC transmit impulse response (N x 1 column vector)
%   gtMatched - Time-reversed RRC matched filter (N x 1 column vector)

if nargin < 5 || isempty(doPlot)
    doPlot = true;
end

f_row = f(:).';
Hf = zeros(size(f_row));

for i = 1:length(f_row)
    if 0 <= abs(f_row(i)) && abs(f_row(i)) < (1-roll_off)/(2*T)
        Hf(i) = T;
    elseif abs(f_row(i)) >= (1-roll_off)/(2*T) && abs(f_row(i)) <= (1+roll_off)/(2*T)
        Hf(i) = T/2 * (1 + cos((pi*T/roll_off) * (abs(f_row(i)) - (1 - roll_off)/(2*T))));
    else
        Hf(i) = 0;
    end
end

ht = ifftshift(ifft(fftshift(Hf)));
Hf = Hf / max(abs(ht));
htNorm = ifftshift(ifft(fftshift(Hf)));
Gf = sqrt(Hf);
gt_row = real(ifftshift(ifft(fftshift(Gf))));

if doPlot
    figure;
    subplot(3, 2, [1, 2]);
    plot(f_row, Hf, 'b', 'LineWidth', 2); hold on;
    plot(f_row, Gf, 'r--', 'LineWidth', 2);
    grid on; xlabel('Frequency [Hz]'); ylabel('Magnitude');
    title('Frequency Domain: RC vs. RRC');
    legend('Raised Cosine H(f)', 'Root Raised Cosine G(f)'); hold off;

    subplot(3, 2, 3);
    plot(t, real(ht), 'g', 'LineWidth', 2);
    grid on; xlabel('Time [s]'); ylabel('Amplitude');
    title('Unnormalized RC h(t)');

    subplot(3, 2, 4);
    plot(t, real(htNorm), 'g', 'LineWidth', 2);
    grid on; xlabel('Time [s]'); ylabel('Amplitude');
    title('Normalized RC h(t)');

    subplot(3, 2, [5, 6]);
    plot(t, gt_row, 'g', 'LineWidth', 2); hold on;
    for shift = -2:2
        if shift ~= 0
            plot(t + shift*T, gt_row, '--', 'LineWidth', 1.2);
        end
    end
    hold off; grid on; xlabel('Time [s]'); ylabel('Amplitude');
    title('Root Raised Cosine Pulse g(t) & Symbol-Spaced Shifts');
    legend('Centered (k=0)', 'Symbol-shifted pulses');
    sgtitle('Half-Nyquist (RRC) Filter Characteristics');
end

gt = gt_row(:);
gtMatched = flipud(gt);
end
