function [toa, cfo] = frame_cfo_acquisition(signal, pilotSymb, K, Tsymb)
%% Joint Time-of-Arrival (ToA) & Carrier Frequency Offset (CFO) Acquisition
% Course: ELEC-H401 (Digital Communications)
%
% Uses differential cross-correlation across lags k = 1..K < length(pilotSymb)
% to locate the known pilot sequence (ToA) robustly under frequency offset and
% estimate CFO [Hz] from the differential phase argument at the peak index.
%
% Inputs:
%   signal    - Downsampled received complex symbol vector
%   pilotSymb - Known complex pilot sequence (length N > K)
%   K         - Differential averaging window length (must satisfy K < N)
%   Tsymb     - Symbol period [s]
%
% Outputs:
%   toa       - Estimated start index of pilot sequence in signal
%   cfo       - Estimated carrier frequency offset [Hz]

y = signal(:);
a = pilotSymb(:);
N = length(a);
K = min(K, N - 1);

D = zeros(length(y), K);
sumK = zeros(length(y), 1);

for k = 1:K
    for n = 1:length(y)-N+1
        l = (k:N-1).';
        sumDk = sum((conj(y(n+l)) .* a(l+1)) .* conj((conj(y(n+l-k)) .* a(l-k+1))), 1);
        D(n, k) = (1 / (N - k)) * sumDk;
    end
    sumK = sumK + abs(D(:, k));
end

[~, toa] = max(sumK);

k_vec = 1:K;
sumDeltaF = sum(angle(D(toa, k_vec)) ./ (2*pi*k_vec*Tsymb));
cfo = -(1/K) * sumDeltaF;
end
