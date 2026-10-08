"""Talc identification against the NIST WebBook talc reference (OH stretching region).

Reference: NIST Chemistry WebBook, NIST Standard Reference Database 69, talc (CAS 14807-96-6),
IR spectrum from Pacific Northwest National Laboratory (public domain), diffuse reflectance.
Only the OH stretching region (3620-3720 cm-1) is used, because the fingerprint region of a
neat-powder diffuse-reflectance spectrum is distorted and not comparable with ATR.

Requires: numpy, pandas.   Run:  python talc_nist_check.py
"""
import numpy as np, pandas as pd

# ---- USER SETTINGS: put this script in the same folder as the data files ----
REF = 'talc_nist_reflectance.txt'   # NIST talc: wavenumber, reflectance (0-1)
DATA = '{}.txt'                     # pan spectra: wavenumber, transmittance (0-1)
PANS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15]
OUT = 'talc_nist_check_results.xlsx'
# ---------------------------------------------------------------------------
W = (3620, 3720)
GRID = np.arange(3625, 3716, 1.0)

def window(x, a):
    m = (x > W[0]) & (x < W[1]); xx, aa = x[m], a[m]
    return xx, aa - np.interp(xx, [xx[0], xx[-1]], [aa[0], aa[-1]])

r = np.loadtxt(REF, comments='#'); rx, ra = window(r[:, 0], np.log10(1 / r[:, 1]))
ref = np.interp(GRID, rx, ra); ref_peak = rx[np.argmax(ra)]

rows = []
for n in PANS:
    d = np.loadtxt(DATA.format(n)); o = np.argsort(d[:, 0])
    px, pa = window(d[o, 0], -np.log10(d[o, 1]))
    best = max((np.corrcoef(np.interp(GRID, px + s, pa), ref)[0, 1], s) for s in np.arange(-6, 6.5, 0.5))
    peak = px[np.argmax(pa)]
    rows.append({'Pan': n, 'Talc OH peak (cm-1)': round(peak, 1), 'NIST talc peak (cm-1)': round(ref_peak, 1),
                 'Shift (cm-1)': round(peak - ref_peak, 1), 'Correlation r (OH region)': round(best[0], 3),
                 'Talc identified (|shift| <= 5)': 'Yes' if abs(peak - ref_peak) <= 5 else 'No'})
T = pd.DataFrame(rows).set_index('Pan'); print(T.to_string()); T.to_excel(OUT)
