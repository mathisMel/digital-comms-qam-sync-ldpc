# Digital Communications — QAM Modulation, Synchronization & Hard-Decision LDPC

A MATLAB baseband physical-layer simulation suite (`ELEC-H401`) modeling **Gray-coded Square-QAM (`4/16/64/256-QAM`) over AWGN**, **carrier/timing synchronization impairments**, **non-data-aided Gardner timing recovery**, **differential pilot ToA/CFO acquisition**, and **rate-1/2 `(128x256)` regular LDPC hard-decision decoding**.

## Modules & Verified Results

1. **Square-QAM & Half-Nyquist Pulse Shaping (`src/step1_qam_ber_awgn.m`)**
   - Constructs a normalized Raised-Cosine / Root-Raised-Cosine (`rolloff = 0.2`, `M = 4`, `N = 100` taps) filter pair via `src/lib/half_nyquist_filter.m`.
   - Evaluates `4/16/64/256-QAM` AWGN BER waterfalls (`Eb/N0 = -5..25 dB`): at `15 dB`, observed BERs are `0`, `0`, `7.67e-4`, and `2.05e-2`.
2. **Constellation Distortion & Impairment BER (`src/step2_sync_impairment_constellations.m`, `src/step3_sync_impairment_ber.m`)**
   - Visualizes constant carrier phase rotation (`phi_0 = pi/8`) and progressive CFO phase drift on 4-QAM.
   - Measures 4-QAM BER degradation (`M = 100`, `Eb/N0 = 0..16 dB`) under CFO ISI (`0..50 ppm`), constant phase shift (`0, pi/8, pi/4` -> `0`, `3.20e-3`, `2.57e-1` @ `10 dB`), and fractional sampling offset (`t_0 = 0..0.30 T`).
3. **Gardner Timing Error Detector (`src/step4_gardner_timing_recovery.m`)**
   - Closed-loop timing lock from `Shift = 0.40 T` across loop gains `K = [0.01, 0.02, 0.03]` (final estimates `[0.3981, 0.4041, 0.4034] T`, final group mean timing error `-3.55e-10..-3.82e-9 s`) and CFO up to `50 ppm`.
4. **Differential Pilot ToA & CFO Acquisition (`src/step5_frame_cfo_acquisition.m`)**
   - Sweeps pilot symbol length `N = [10, 20, 40]` (`K = 8`) and averaging window `K = [1, 8, 16]` (`N = 20`) over `20` trials (`10 ppm` CFO, `600 MHz`): at `16 dB`, all lengths achieve `0.00` mean ToA error with CFO stds `[10.76, 4.66, 1.78] ppm`.
5. **End-to-End Synchronized 16-QAM Receiver (`src/step6_full_sync_receiver.m`)**
   - Chains `Shift = 0.17 T` Gardner lock, pilot acquisition (`toa = 930`, `est_CFO = 30.11 kHz` vs. `30.00 kHz` at `50 ppm`), and unwrapped linear phase interpolation — yielding `0` bit errors across `0, 2, 10 ppm`.
6. **Hard-Decision LDPC Coding (`src/step7a_ldpc_small_matrix_demo.m`, `src/step7b_ldpc_hard_ber.m`)**
   - Verifies systematic Gauss-Jordan encoding on a `(10, 5)` matrix (`0` mismatches) and simulates `(128x256)` rate-1/2 LDPC hard-decision Tanner-graph decoding (`1, 2, 5, 10` iterations) from `-4..10 dB` (`BER = 0.398 -> 3.26e-4 -> 0` @ `8 dB`).
   - *Scope note:* Soft-decision (LLR / belief propagation) LDPC decoding (`Report Section 3.3`) is omitted because only hard-decision decoders were present in the course source archive.

## Repository Structure

```text
digital-comms-qam-sync-ldpc/
├── src/
│   ├── step1_qam_ber_awgn.m
│   ├── step2_sync_impairment_constellations.m
│   ├── step3_sync_impairment_ber.m
│   ├── step4_gardner_timing_recovery.m
│   ├── step5_frame_cfo_acquisition.m
│   ├── step6_full_sync_receiver.m
│   ├── step7a_ldpc_small_matrix_demo.m
│   ├── step7b_ldpc_hard_ber.m
│   └── lib/                     # mapping, demapping, RRC filter, AWGN, Gardner, ToA/CFO, LDPC
├── docs/figures/                # Verified exported MATLAB plots
├── .gitignore
├── LICENSE
├── THIRD_PARTY_NOTICES.md       # Attribution for Bagawan Nugroho LDPC routines
└── README.md
```

## MATLAB Requirements & Quickstart

- **MATLAB + Communications Toolbox** (`bi2de` in `src/lib/mapping.m`; verified on MATLAB R2026b).
```matlab
run('src/step1_qam_ber_awgn.m')
run('src/step3_sync_impairment_ber.m')
run('src/step5_frame_cfo_acquisition.m')
run('src/step7b_ldpc_hard_ber.m')
```

## Key Figures

![Square-QAM BER over AWGN](docs/figures/step1_qam_fig03.png)
![4-QAM BER under Constant Carrier Phase Shift](docs/figures/step3_imp_fig03.png)
![Gardner Timing Error Convergence](docs/figures/step4_gardner_fig03.png)
![CFO Estimation Std vs. Pilot Length](docs/figures/step5_acq_fig03.png)
![Rate-1/2 Hard-Decision LDPC BER](docs/figures/step7b_ldpc_fig01.png)
