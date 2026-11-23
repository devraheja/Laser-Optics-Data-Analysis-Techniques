# Nd:YAG results

## Pump diode (fit to I > 200 mA)
- slope s = 1.0661 ± 0.0052 W/A, intercept b = -213.02 ± 2.17 mW
- threshold I_th = -b/s = 199.8 mA
  - error ignoring the s-b correlation: ± 2.26 mA
  - error with full covariance: ± 1.14 mA (fitting I_th directly gives ± 1.14 mA)
  - residual bootstrap 68% interval: [198.8, 200.8] mA
- fit-window sensitivity:

| lowest current (mA) | points | I_th (mA) |
|---|---|---|
| 225 | 9 | 199.8 ± 1.1 |
| 250 | 8 | 199.9 ± 1.5 |
| 300 | 7 | 201.1 ± 2.0 |
| 350 | 6 | 199.7 ± 3.0 |
| 400 | 5 | 194.4 ± 3.0 |

## Laser characteristic
- linear fit: slope efficiency η = 0.120 ± 0.008, P_th = 116.5 ± 11.9 mW (fit 11.9, pump calibration 0.8)
- quadratic fit: P_th = 77.0 mW
- data bracket: no lasing at 80.2 mW, lasing at 106.8 mW
- residual signs of the linear fit: + - - - - + +
