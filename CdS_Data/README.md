# Raw experimental data

The nine files in this directory are the only inputs to the pipeline.
Everything else in the repository is derived from them.

All files are read by `GH_Data_Processing.m` and are never modified. Do not
re-save them from a spreadsheet application: the CSV files use `;` as the
field separator and `.` as the decimal mark, and an Italian-locale round
trip through Excel will silently swap the two.

## Sulfur speciation maps

Quantitative sulfur phase maps obtained by linear combination fitting (LCF)
of full-spectral S K-edge µ-XANES imaging stacks, acquired at ESRF-ID21
(experiment HG-227) on 15 µm thick sections of hexagonal CdS oil paint
mock-ups. Values are relative concentrations in [0, 1].

| File | Sample | Map ID | Exposure | Grid (rows × cols) | Columns used |
|---|---|---|---|---|---|
| `Sample_t1.csv` | t1 | C2 | 24 h | 194 × 14 | 3–4 |
| `Sample_t4.csv` | t4 | C1 | 80 h | 284 × 17 | 3–4 |
| `Sample_t7.csv` | t7 | D1 | 167 h | 156 × 15 | 4–5 |
| `Sample_t9.csv` | t9 | B2 | 340 h | 187 × 20 | 3–4 |
| `Sample_t10.csv` | t10 | B1 | 680 h | 234 × 15 | 3–4 |

The *Map ID* column is the internal identifier of the XANES map and appears
throughout the MATLAB code as the variable suffix (`CdSC2_crop`,
`Adim_D1_mean`, `fitted_value_B1`, and so on). The mapping between sample
name and map ID is fixed and is the single most important thing to keep in
mind when reading the code:

```
t1 -> C2      t4 -> C1      t7 -> D1      t9 -> B2      t10 -> B1
```

**The column layout is not uniform across files.** Sample t7 carries an
extra leading intensity column, so its CdSO4 and CdS fits sit in columns
4–5 rather than 3–4. This is why `GH_Data_Processing.m` reads
`MD(:,4:5)` for t7 and `M(:,3:4)` for the others. Sample t4 also has a
trailing intensity column, which is simply ignored.

Common structure of every file:

- `row`, `column` — pixel indices of the map, row-major
- `CdSO4-fit` — relative CdSO4 concentration from the LCF
- `CdS-fit` (`hex-CdS_fit` in t9) — relative hexagonal CdS concentration

Only the CdS column feeds the model; the CdSO4 column is read into the
`CdSO4*` variables and retained for reference. Pixel size is 0.25 µm along
the depth (`dx`) and 1 µm laterally (`dy`), both hard-coded in
`GH_Data_Processing.m`.

## Lamp emission spectra

| File | Sample |
|---|---|
| `Lamp_Spectrum_t4.txt` | t4 |
| `Lamp_Spectrum_t7.txt` | t7 |
| `Lamp_Spectrum_t10.txt` | t10 |

Tab-separated, two header lines (`Wavelength / Media / StDev`, then units),
about 1625 data rows. Columns: wavelength [nm], mean irradiance, standard
deviation. Only the first two columns are used.

Spectra exist for three of the five samples. The pipeline truncates each at
the CdS band-gap wavelength (512 nm), averages the three, and normalizes by
the maximum to obtain the dimensionless illumination profile
`I_mean_adim` on the dimensionless wavelength grid `lambda_adim`. The
manuscript notes that irradiance varied between experiments and that the
model is driven by a single representative average (ca. 1.33 × 10^5
µW/cm²); this averaging step is where that simplification enters the code.

Note that early rows contain small negative irradiance values (detector
baseline). They survive the average and are clipped later by `subplus`
inside the solvers, not here.

## UV-Vis reflectance

`Reflectance_UV_Vis.txt` — tab-separated, no header, 651 rows. Column 1 is
wavelength [nm], column 2 is reflectance [%].

Wavelengths run **descending**, from 850 nm down to 200 nm, which is why
`GH_Data_Processing.m` reverses the matrix (`R(end:-1:1,:)`) before use.
Reflectance is divided by 100 and made non-negative to give `Refl_un`,
from which the solvers derive the CdS molar absorptivity through

```
eps_c_values = -log(Refl_un) / (c_c_ref * L)
```

with `c_c_ref = 4.82/144.46` mol/cm³ and `L` the layer thickness. The
CdSO4 absorptivity is obtained as `eps_c_values / nu`, so the reflectance
of a single reference material drives both species, with `nu` calibrated.

## Provenance and citation

These files are a working subset. The complete deposited datasets are:

- S K-edge XANES data — ESRF-ICAT: <https://doi.esrf.fr/10.15151/ESRF-DC-2491544156>
- External reflection FTIR and VIS spectroscopy — Zenodo: <https://doi.org/10.5281/zenodo.21652776>

If you use these data, cite the manuscript and the deposited datasets, not
only this repository. See the root `README.md` for the citation block.
