function [bit_rx] = demapping(symb_rx, Nbps, modulation)
%% Hard-Decision Gray PAM / Square-QAM Demapper
% Course: ELEC-H401 (Digital Communications)
%
% Inverse of mapping.m for hard-decision detection of unit-variance PAM or
% square-QAM constellations. Notice bit extraction writes MSB-first directly
% into column 1..NbpsI/Q without extra fliplr reversals.
%
% Inputs:
%   symb_rx    - Column vector of received equalized symbols (variance ~1)
%   Nbps       - Bits per symbol (2,4,6,8 for square QAM)
%   modulation - 'pam' or 'qam'
%
% Outputs:
%   bit_rx     - Recovered binary column vector of length Nsymb * Nbps

symb_rx = symb_rx(:);
Nsymb = size(symb_rx, 1);

switch modulation
    case 'pam'
        sigma = sqrt(sum(((0:2^Nbps-1) - (2^Nbps-1)/2).^2) / 2^Nbps);
        int_rx = sigma * symb_rx + (2^Nbps-1)/2;

        int_det = round(int_rx);
        int_det(int_det < 0) = 0;
        int_det(int_det > 2^Nbps-1) = 2^Nbps-1;

        k = size(int_det, 2);
        weights = 2.^(k-1:-1:0);
        mapp_rx = int_det * weights.';

        bit_rx2(:, 1) = mapp_rx(:, 1);
        for ii = 2:Nbps
            bit_rx2(:, ii) = xor(mapp_rx(:, ii-1), mapp_rx(:, ii));
        end
        bit_rx = reshape(bit_rx2', Nsymb*Nbps, 1);

    case 'qam'
        % Real part (I)
        NbpsI = Nbps/2;
        symb_rxI = real(symb_rx);

        sigmaI = sqrt(sum(((0:2^NbpsI-1) - (2^NbpsI-1)/2).^2) / 2^NbpsI);
        int_rxI = sigmaI * sqrt(2) * symb_rxI + (2^NbpsI-1)/2;

        int_detI = round(int_rxI);
        int_detI(int_detI < 0) = 0;
        int_detI(int_detI > 2^NbpsI-1) = 2^NbpsI-1;

        mapp_rxI = zeros(Nsymb, NbpsI);
        for i = 1:Nsymb
            q = int_detI(i);
            for k = 1:NbpsI
                mapp_rxI(i, NbpsI-k+1) = mod(q, 2);
                q = floor(q/2);
            end
        end

        bit_rx2I(:, 1) = mapp_rxI(:, 1);
        for ii = 2:NbpsI
            bit_rx2I(:, ii) = xor(mapp_rxI(:, ii-1), mapp_rxI(:, ii));
        end

        % Imaginary part (Q)
        NbpsQ = Nbps/2;
        symb_rxQ = imag(symb_rx);

        sigmaQ = sqrt(sum(((0:2^NbpsQ-1) - (2^NbpsQ-1)/2).^2) / 2^NbpsQ);
        int_rxQ = sigmaQ * sqrt(2) * symb_rxQ + (2^NbpsQ-1)/2;

        int_detQ = round(int_rxQ);
        int_detQ(int_detQ < 0) = 0;
        int_detQ(int_detQ > 2^NbpsQ-1) = 2^NbpsQ-1;

        mapp_rxQ = zeros(Nsymb, NbpsQ);
        for i = 1:Nsymb
            q = int_detQ(i);
            for k = 1:NbpsQ
                mapp_rxQ(i, NbpsQ-k+1) = mod(q, 2);
                q = floor(q/2);
            end
        end

        bit_rx2Q(:, 1) = mapp_rxQ(:, 1);
        for ii = 2:NbpsQ
            bit_rx2Q(:, ii) = xor(mapp_rxQ(:, ii-1), mapp_rxQ(:, ii));
        end

        bit_rx = reshape([bit_rx2I, bit_rx2Q]', Nsymb*Nbps, 1);
end
end
