function vNodes = ldpc_hard_decode(r, H, maxIter)
%% Hard-Decision LDPC Bit-Flipping / Message-Passing Decoder
% Course: ELEC-H401 (Digital Communications)
%
% Implements iterative hard-decision Tanner-graph decoding with majority voting
% until parity-check syndrome H*v' == 0 (mod 2) or maxIter is reached.
%
% Inputs:
%   r       - Received binary codeword row vector (1 x N)
%   H       - Parity-check matrix (K x N)
%   maxIter - Maximum decoding iterations
%
% Outputs:
%   vNodes  - Decoded binary codeword row vector (1 x N)

[K, N] = size(H);
vNodes = r;
vToC = zeros(K, N);
syndrome = vNodes * H';
iter = 0;

while iter < maxIter && norm(mod(syndrome, 2)) ~= 0
    % Step 0: Variable nodes to check nodes
    for c = 1:N
        index = find(H(:, c));
        vToC(index, c) = vNodes(c);
    end

    % Step 1: Check-node parity updates
    cToV = zeros(K, N);
    for k = 1:K
        index = find(H(k, :));
        for i = 1:length(index)
            newIndex = index;
            newIndex(i) = [];
            cToV(k, index(i)) = mod(sum(vToC(k, newIndex)), 2);
        end
    end

    % Step 2: Variable-node majority vote
    for c = 1:N
        index = find(H(:, c));
        vote = (sum(cToV(index, c)) + vNodes(c)) / (length(index) + 1);
        vNodes(c) = double(vote > 0.5);
    end

    syndrome = vNodes * H.';
    iter = iter + 1;
end
end
