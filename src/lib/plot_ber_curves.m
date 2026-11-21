function plot_ber_curves(SNR, Y, legends, yLabelStr, titleStr, useLog)
%% Standardized Multi-Curve Performance Plotter
% Course: ELEC-H401 (Digital Communications)
%
% Inputs:
%   SNR       - Horizontal E_b/N_0 vector [dB]
%   Y         - Performance matrix (numCurves x numel(SNR))
%   legends   - Cell array of legend labels
%   yLabelStr - Optional y-axis label (default 'BER')
%   titleStr  - Optional plot title (default 'SNR vs. BER')
%   useLog    - Optional boolean for log y-axis (default true)

if nargin < 4 || isempty(yLabelStr), yLabelStr = 'BER'; end
if nargin < 5 || isempty(titleStr), titleStr = 'SNR vs. BER'; end
if nargin < 6 || isempty(useLog), useLog = true; end

figure;
hold on;
for i = 1:size(Y, 1)
    plot(SNR, Y(i, :), '-o', 'LineWidth', 1.2);
end
if useLog
    set(gca, 'YScale', 'log');
end
grid on;
xlabel('E_b/N_0 [dB]');
ylabel(yLabelStr);
title(titleStr);
legend(legends, 'Location', 'best');
hold off;
end
