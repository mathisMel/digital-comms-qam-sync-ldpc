function [symbolsCorrected, epsilon] = gardner_ted(symbols, K, M)
%% Non-Data-Aided Gardner Timing Error Detector & Closed-Loop Recovery
% Course: ELEC-H401 (Digital Communications)
%
% Tracks fractional symbol timing offset epsilon(n) using mid-symbol
% linear interpolation and updates the sampling instant iteratively.
%
% Inputs:
%   symbols - Oversampled matched-filter output column vector (rate M/T)
%   K       - Loop gain step size
%   M       - Oversampling ratio (samples per symbol)
%
% Outputs:
%   symbolsCorrected - Timing-recovered symbol-rate column vector
%   epsilon          - Estimated fractional symbol timing offset per symbol

epsilon = zeros(floor(length(symbols)/M), 1);
symbolsCorrected = zeros(size(epsilon));
indiceLowStock = zeros(size(epsilon));
indiceHighStock = zeros(size(epsilon));

symbolsCorrected(1) = symbols(1);
epsilon(1) = 0;

for n = 1:length(epsilon)-1
    indiceLow = floor(((n + epsilon(n))*M) + 1);
    indiceLowStock(n+1) = indiceLow;
    indiceHigh = ceil(((n + epsilon(n))*M) + 1);
    indiceHighStock(n+1) = indiceHigh;
    symbolsCorrected(n+1) = (symbols(indiceHigh) + symbols(indiceLow)) / 2;

    middleIndiceHigh = ceil((indiceHighStock(n) + ((n + epsilon(n))*M) + 1) / 2);
    middleIndiceLow  = floor((indiceLowStock(n) + ((n + epsilon(n))*M) + 1) / 2);
    middlesymbol = (symbols(middleIndiceHigh) + symbols(middleIndiceLow)) / 2;

    epsilon(n+1) = epsilon(n) - K * real(middlesymbol * (conj(symbolsCorrected(n+1)) - conj(symbolsCorrected(n))));
end
end
