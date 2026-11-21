%% Digital Comms — Step 7b: Rate-1/2 (128x256) LDPC Hard-Decision BER over AWGN
% Course: ELEC-H401 (Digital Communications)
%
% Evaluates rate-1/2 regular LDPC coding with 4-QAM modulation and iterative
% hard-decision Tanner-graph decoding (1, 2, 5, 10 iterations) across Eb/N0.
% Uses Bagawan Nugroho's generate_ldpc / encode_ldpc helpers (see THIRD_PARTY_NOTICES.md).
% Note: Soft-decision (LLR) LDPC from Report Section 3.3 is not included in this archive.

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
RollOff = 0.2;
Rs = 5e6;
M = 4;
N = 100;
Hrow = 128;
Hcolumn = 256;
SNR = -4:2:10;              % Eb/N0 [dB]
Iterations = [1, 2, 5, 10];
Nexp = 3;                   % Blocks of 2048 info bits averaged per point

T = 1/Rs;
Fs = Rs*M;
df = (1/N)*Fs;
f = -Fs/2:df:((Fs/2)-df);
Ts = 1/Fs;
t = (-N/2:(N/2)-1)*Ts;
Nb = 2048;

BER = zeros(numel(Iterations), numel(SNR));
decodedData = zeros(Nb/Hrow, Hcolumn);

%% Filter, Parity-Check Matrix & Encoded Frame
[gt, gtMatched] = half_nyquist_filter(T, RollOff, f, t, false);
H0 = generate_ldpc(Hrow, Hcolumn, 0, 1, 3);
inputBitVec = randi([0 1], [Nb, 1]);
inputBitMatrixBlock = reshape(inputBitVec, Hrow, Nb/Hrow);

[parityBitsMatrixBlock, H] = encode_ldpc(inputBitMatrixBlock, H0, 0);
codewordBlocks = [parityBitsMatrixBlock; inputBitMatrixBlock];
parityBitVec = reshape(codewordBlocks, [], 1);

symbTx = mapping(parityBitVec, Nbps, modulation);
symbTxUp = upsample(symbTx, M);
symbTxUpFiltered = conv(symbTxUp, gt);
Eb = compute_eb(symbTx, Nbps, T);

%% AWGN BER Sweep Across Decoder Iterations
for l = 1:Nexp
    for i = 1:length(Iterations)
        for j = 1:length(SNR)
            symbTxCorrupted = add_awgn(Eb, SNR(j), Rs, symbTxUpFiltered);
            symbRxCorruptedFiltered = conv(symbTxCorrupted, gtMatched);
            symbRxCorruptedCropped = symbRxCorruptedFiltered(N:end-(N-1));
            symbRxCorruptedFilteredDown = downsample(symbRxCorruptedCropped, M);
            bitRxCorrupted = demapping(symbRxCorruptedFilteredDown, Nbps, modulation);
            outputBitMatrixBlock = reshape(bitRxCorrupted, Hcolumn, Nb/Hrow).';

            for k = 1:Nb/Hrow
                decodedData(k, :) = ldpc_hard_decode(outputBitMatrixBlock(k, :), H, Iterations(i));
            end

            decodedDataForError = decodedData(:, (Hcolumn-Hrow+1):end).';
            error = numel(inputBitMatrixBlock) - nnz(inputBitMatrixBlock == decodedDataForError);
            BER(i, j) = BER(i, j) + error / Nb;
        end
    end
end
BER = BER ./ Nexp;

%% Plot BER Curves
plot_ber_curves(SNR, BER, {'1 iteration', '2 iterations', '5 iterations', '10 iterations'}, ...
    'BER', 'Rate-1/2 (128x256) Hard-Decision LDPC BER (4-QAM)', true);
ylim([1e-4 1]);

fprintf('Step 7b LDPC BER (10 iter) across SNR=-4..10 dB: %s\n', sprintf('%.3e ', BER(4,:)));
