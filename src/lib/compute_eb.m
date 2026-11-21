function [E_b] = compute_eb(symb_tx, Nbps, T)
%% Energy-per-Bit Calculator for Passband Equivalent Signals
% Course: ELEC-H401 (Digital Communications)
%
% Computes passband energy per bit E_b = (E[|s_bb|^2] / 2) * (T / Nbps).
%
% Inputs:
%   symb_tx - Complex baseband symbol vector
%   Nbps    - Bits per symbol
%   T       - Symbol duration [s]
%
% Outputs:
%   E_b     - Passband energy per bit [J/bit]

var_bb_signal = (1 / numel(symb_tx)) * sum(abs(symb_tx).^2);
var_signal = var_bb_signal / 2; % Factor 1/2 from passband real projection
T_b = T / Nbps;
E_b = var_signal * T_b;
end
