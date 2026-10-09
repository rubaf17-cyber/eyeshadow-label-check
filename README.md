# Eyeshadow Label Check
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23250764.svg)](https://doi.org/10.5281/zenodo.23250764)

Open-source ATR-FTIR tool for screening eyeshadow pans for eight common base ingredients and checking them against label claims.

**Live tool:** https://rubaf17-cyber.github.io/eyeshadow-label-check/

## What it does
- Reads an ATR-FTIR spectrum (wavenumber vs. %T, T or absorbance)
- Screens eight base ingredients (talc, mica, kaolin, PDMS, ester emollients, metal stearate, MgCO₃, CaCO₃) by their diagnostic bands, using transparent rule-based criteria (Savitzky–Golay second-derivative signal-to-noise)
- Identifies the stearate metal from the carboxylate band position (Zn ≈ 1540 cm⁻¹; Mg ≈ 1575 cm⁻¹; Ca 1577 + 1540 cm⁻¹)
- Confirms talc by comparing the OH stretching band (3620–3720 cm⁻¹) with the NIST talc reference (peak position and correlation)
- The user compares the result with the label

## Diagnostic bands and references
| Component | Band (cm⁻¹) | Reference |
|---|---|---|
| Talc | 3674 (OH), 667 | NIST Chemistry WebBook, CAS 14807-96-6 (PNNL) |
| Mica | 3626 (ratio to 3674) | NIST / Coblentz No. 4666 |
| Kaolin | 3689 | Madejová & Komadel (2001) |
| PDMS (dimethicone) | 1260, 800 | NIST Chemistry WebBook, CAS 63148-62-9 (PNNL) |
| Ester emollients | 1740 | McEwan et al. (2025) |
| Zinc stearate | 1535–1540 | Hermans & Helwig (2020); Larkin & Jackson (2024); NIST / Coblentz No. 1069 |
| Magnesium stearate | 1573–1578 | SDBS No. 12704 (Nujol and KBr disk) |
| Magnesium carbonate | 1446, 885, 748 | Weir & Lippincott (1961) |
| Calcium carbonate | 1420, 872, 712 | Weir & Lippincott (1961); NIST / Coblentz Nos. 4659, 4660 |

## Files
| File | Description |
|---|---|
| `index.html` | Web tool (runs in any browser, no installation) |
| `eyeshadow_ftir_screening.py` | Rule-based screening, PCA and HCA (Python) |
| `talc_nist_check.py` | Talc identification against the NIST talc reference (Python) |
| `talc_nist_reflectance.txt` | NIST talc reference spectrum (wavenumber, reflectance); used by the tool for the OH-region match |
| `pdms_nist_k.txt` | NIST PDMS reference (wavenumber, imaginary refractive index k) |
| `reference_bands.csv` | All diagnostic bands used by the tool, with assignment, source and link |

## Reference data
| Material | Source | In this repository |
|---|---|---|
| Talc, CAS 14807-96-6 | NIST WebBook, PNNL (public domain): https://webbook.nist.gov/cgi/cbook.cgi?ID=14807-96-6&Units=SI&cIR=on | `talc_nist_reflectance.txt` |
| PDMS, CAS 63148-62-9 | NIST WebBook, PNNL (public domain): https://webbook.nist.gov/cgi/cbook.cgi?ID=63148-62-9&Units=SI&cIR=on | `pdms_nist_k.txt` |
| Mica, Coblentz No. 4666 | NIST WebBook, Coblentz Society: https://webbook.nist.gov/cgi/cbook.cgi?ID=B6004666&Units=SI&cIR=on | link only (Coblentz copyright) |
| Calcite, Coblentz No. 4659 | NIST WebBook, Coblentz Society: https://webbook.nist.gov/cgi/cbook.cgi?ID=B6004659&Units=SI&cIR=on | link only (Coblentz copyright) |
| CaCO₃ precipitated, CAS 471-34-1, Coblentz No. 4660 | NIST WebBook, Coblentz Society: https://webbook.nist.gov/cgi/cbook.cgi?ID=C471341&Units=SI&cIR=on | link only (Coblentz copyright) |
| Kaolin | Madejová & Komadel (2001), ATR band table: https://doi.org/10.1346/CCMN.2001.0490508 | band positions in `reference_bands.csv` |
| Ester emollients | McEwan et al. (2025), *Forensic Chem.* 45, 100677, Table 3 | band positions in `reference_bands.csv` |
| Zinc stearate | Hermans & Helwig (2020): https://doi.org/10.1177/0003702820935183 ; Larkin & Jackson (2024): https://doi.org/10.1177/27551857241253834 | band positions in `reference_bands.csv` |
| Magnesium stearate | SDBS No. 12704, Nujol and KBr disk (AIST): https://sdbs.db.aist.go.jp/CompoundLanding.aspx?sdbsno=12704 | band positions only (SDBS data may not be redistributed) |
| MgCO₃ and CaCO₃ | Weir & Lippincott (1961): https://doi.org/10.6028/jres.065A.021 | band positions in `reference_bands.csv` |

The PNNL files were parsed from the JCAMP-DX files downloaded from the NIST WebBook (header lists the owner as public domain).

## Requirements (Python scripts)
numpy, scipy, pandas, matplotlib, openpyxl

## Limitations
- Detection thresholds were set on 14 pans measured on one diamond-ATR instrument (one spectrum per pan); validate on your own instrument.
- "Not detected" means below the detection rule, not proven absent.
- The NIST talc spectrum was measured by diffuse reflectance of neat powder; only its OH stretching region is comparable with ATR and is used here.
- Synthetic mica has no OH band and is not seen in the OH region.
- Pigments are not tested.

## References
1. NIST Chemistry WebBook, NIST Standard Reference Database Number 69: talc (CAS 14807-96-6) and polydimethylsiloxane (CAS 63148-62-9), data from Pacific Northwest National Laboratory. https://webbook.nist.gov
2. NIST / Coblentz Society IR spectra collection: Nos. 4666 (mica), 4659 and 4660 (calcite), 1069 (zinc stearate).
3. Madejová, J.; Komadel, P. (2001). Baseline studies of the Clay Minerals Society source clays: infrared methods. *Clays and Clay Minerals* 49, 410–432. https://doi.org/10.1346/CCMN.2001.0490508
4. McEwan, et al. (2025). *Forensic Chemistry* 45, 100677.
5. SDBS, Spectral Database for Organic Compounds, AIST, Japan. No. 12704 (magnesium stearate). https://sdbs.db.aist.go.jp
6. Hermans, J.; Helwig, K. (2020). The identification of multiple crystalline zinc soap structures using infrared spectroscopy. *Applied Spectroscopy* 74(12). https://doi.org/10.1177/0003702820935183
7. Larkin, P. J.; Jackson, A. (2024). Interpretation of the infrared spectra of metal-stearate salts. *Applied Spectroscopy Practica*. https://doi.org/10.1177/27551857241253834
8. Weir, C. E.; Lippincott, E. R. (1961). Infrared studies of aragonite, calcite, and vaterite type structures in the borates, carbonates, and nitrates. *J. Res. Natl. Bur. Stand.* 65A, 173–183. https://doi.org/10.6028/jres.065A.021

## Version
v0.4: upload-only interface (example pans removed); the spectrum is analysed in the browser and not uploaded anywhere.
v0.3: label input removed; kaolin band moved to 3689 cm⁻¹; stearate metal type added; reference list added.

## Citation
Abbas, R. F. (2026). Eyeshadow Label Check: an open-source ATR-FTIR screening tool for eyeshadow base ingredients (v0.4.1). Zenodo. https://doi.org/10.5281/zenodo.23250764

## License
MIT
