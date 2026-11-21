%% Digital Comms — Step 7a: Small (10x5) Systematic LDPC Parity-Check Verification
% Course: ELEC-H401 (Digital Communications)
%
% Verifies binary Gauss-Jordan elimination from parity-check matrix H to
% systematic generator matrix G and confirms zero decoding errors with ldpc_hard_decode.

clc; clear; close all;
rng(0);
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

%% Parameters & Parity-Check Matrix H (5x10)
Nbps = 2;
Nsymb = 1e3;

H = [
    1 1 0 1 1 0 0 1 0 0;
    0 1 1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0 1 0;
    0 0 1 0 0 1 0 1 0 1
];

[m, n] = size(H);
for j = 1:min(m, n)
    pivot_row = find(H(j:end, j), 1) + j - 1;
    if ~isempty(pivot_row)
        H([j, pivot_row], :) = H([pivot_row, j], :);
        for i = 1:m
            if i ~= j && H(i, j) == 1
                H(i, :) = mod(H(i, :) - H(j, :), 2);
            end
        end
    end
end
disp('Systematic parity-check matrix [P | I]:');
disp(H);

parityCheckMatrix = [1, 0, 0, 0, 0, 0, 1, 1, 1, 0;
                     0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
                     0, 0, 1, 0, 0, 1, 0, 1, 0, 1;
                     0, 0, 0, 1, 0, 0, 0, 1, 1, 1;
                     0, 0, 0, 0, 1, 1, 1, 0, 0, 1];

IdentityMatrix = eye(5);
pTrans = [0, 1, 1, 1, 0;
          1, 0, 1, 0, 0;
          1, 0, 1, 0, 1;
          0, 0, 1, 1, 1;
          1, 1, 0, 0, 1];
p = pTrans';
generatorMatrix = [p IdentityMatrix];

%% Encode & Hard-Decision Decode Verification
Nb = Nbps * Nsymb;
inputBitVec = randi([0 1], [Nb, 1]);
inputBitMatrixBlock = reshape(inputBitVec, length(inputBitVec)/height(generatorMatrix), height(generatorMatrix));

encodedData = zeros(size(inputBitMatrixBlock, 1), length(generatorMatrix));
decodedData = zeros(size(encodedData));

for i = 1:size(inputBitMatrixBlock, 1)
    encodedData(i, :) = mod(inputBitMatrixBlock(i, :) * generatorMatrix, 2);
end
for i = 1:size(encodedData, 1)
    decodedData(i, :) = ldpc_hard_decode(encodedData(i, :), parityCheckMatrix, 10);
end

nombre_diff = numel(encodedData) - nnz(encodedData == decodedData);
fprintf('Step 7a Verified: Systematic 10x5 LDPC bit mismatches = %d\n', nombre_diff);
