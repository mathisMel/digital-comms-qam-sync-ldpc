function [symb_corrupted, n] = add_awgn(E_b, SNR_db, bw_ref, symb_tx_up_filtered)
%% Complex Baseband Additive White Gaussian Noise (AWGN) Channel
% Course: ELEC-H401 (Digital Communications)
%
% Adds circularly symmetric complex Gaussian noise calibrated to E_b/N_0.
% Note: Because half_nyquist_filter normalizes sum(|gt|^2) = 1 at the matched
% filter output, passing bw_ref = Rs yields exact E_b/N_0 after matched filtering
% regardless of oversampling factor M.
%
% Inputs:
%   E_b                  - Passband energy per bit [J/bit] (from compute_eb)
%   SNR_db               - Target E_b/N_0 [dB]
%   bw_ref               - Noise reference bandwidth [Hz] (use Rs for exact Eb/N0)
%   symb_tx_up_filtered  - Pulse-shaped baseband signal vector
%
% Outputs:
%   symb_corrupted       - Noisy received baseband signal vector
%   n                    - Complex AWGN sample vector

N_0 = E_b / (10^(SNR_db/10));
n_1 = sqrt(N_0 * bw_ref * 2) * randn(size(symb_tx_up_filtered));
n_2 = sqrt(N_0 * bw_ref * 2) * randn(size(symb_tx_up_filtered));
n = sqrt(1/2) * (n_1 + 1j*n_2);
symb_corrupted = symb_tx_up_filtered + n;
end
