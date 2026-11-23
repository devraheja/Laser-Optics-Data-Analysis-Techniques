# Data

All numbers were measured in the Advanced Physics Lab course at FAU Erlangen-Nürnberg and are transcribed from the tables of my lab reports. The raw oscilloscope traces, knife-edge scans and camera images were not kept, so this folder holds the processed tables only.

| File | Experiment | Content |
|---|---|---|
| `hene_knife_edge.csv` | HeNe laser | Parameters of the erf fits to the knife-edge scans: background, amplitude, centre and 1/e² radius with fit errors, at three distances for two resonator lengths |
| `hene_clip_level.csv` | HeNe laser | Blade positions at 90 % and 10 % transmitted power (clip level ε = 0.1). Position reading error ±20 µm |
| `hene_stability.csv` | HeNe laser | Photodiode voltage vs. mirror separation near the stability limit. Background 0.25 V |
| `yag_pump_diode.csv` | Nd:YAG laser | Pump diode optical power vs. current |
| `yag_laser_output.csv` | Nd:YAG laser | 1064 nm output behind the OC100 coupler vs. diode current. No output for I ≤ 275 mA |
| `michelson_calibration.csv` | Michelson / FTS | FFT peak index of the 528 nm and 653 nm interference filters |
| `michelson_slit_peaks.csv` | Michelson / FTS | FFT peak index of eight prism-selected bands of a halogen lamp |
| `michelson_double_slit.csv` | Michelson / FTS | FFT peak pairs for five double-slit configurations (empty: not resolved) |

Also used from the Michelson report: white-light spectral FWHM 26 bins, fringe-envelope FWHM 41 px, spectral barycenter ≈ 735 nm, single-band width 4 to 5 bins, camera rows of 1280 px.
