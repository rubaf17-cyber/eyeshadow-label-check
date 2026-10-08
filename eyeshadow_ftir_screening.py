"""
Rule-based ATR-FTIR screening of eyeshadow pans + exploratory PCA/HCA.
Input : one two-column txt file per pan (wavenumber, transmittance fraction).
Output: band S/N table, rule-based component calls, PCA/HCA figure.

Band intensity = -(2nd derivative of absorbance) at the band (Savitzky-Golay,
window 11, poly 3), divided by the noise of the 2nd derivative in the
featureless 2450-2650 cm-1 region (S/N). Thresholds were set empirically on
this data set and must be validated with reference standards before use
on new samples.
"""
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "{}.txt"   # put the script in the same folder as 1.txt, 2.txt ...
PANS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15]
PALETTE = {p: ("Palette A (1-8)" if p <= 8 else "Haya Queen (9-15)") for p in PANS}

def load(n):
    d = np.loadtxt(DATA.format(n)); o = np.argsort(d[:, 0])
    return d[o, 0], -np.log10(d[o, 1])

# ---------------- 1. band features ----------------
BANDS = {"3674 talc OH": 3674, "3689 kaolin OH": 3689, "3626 Al-OH": 3626, "2917 CH2": 2917,
         "1740 C=O ester": 1740, "1540 COO- Zn": 1540, "1576 COO- Mg": 1576, "1446 CO3 nu3": 1446,
         "1260 Si-CH3": 1260, "872 CO3 nu2": 872, "747 CO3 nu4/Si-O-Al": 748,
         "712 CaCO3 nu4": 712, "667 talc Mg-OH": 667}

def local_prom(x, A, c, w, l, r):
    m = (x > c - w) & (x < c + w); i = np.argmax(np.where(m, A, -9))
    al, ar = np.interp(l, x, A), np.interp(r, x, A)
    return A[i] - (al + (ar - al) * (x[i] - l) / (r - l))

rows = []
for n in PANS:
    x, A = load(n)
    d2 = savgol_filter(A, 11, 3, deriv=2)
    noise = np.std(d2[(x > 2450) & (x < 2650)])
    f = {k: -d2[(x > c - (4 if c == 3689 else 6)) & (x < c + (4 if c == 3689 else 6))].min() / noise
         for k, c in BANDS.items()}
    f["mica ratio 3626/3674"] = (local_prom(x, A, 3626, 6, 3605, 3650) /
                                 local_prom(x, A, 3674, 3, 3660, 3690))
    rows.append({"Pan": n, "Palette": PALETTE[n], **f})
F = pd.DataFrame(rows).set_index("Pan")

# ---------------- 2. rules ----------------
def lvl(v, hi, lo, tr=None):
    if v >= hi: return "Detected"
    if v >= lo: return "Minor" if tr is None else "Minor"
    if tr is not None and v >= tr: return "Trace"
    return "Not detected"

calls = []
for n, r in F.iterrows():
    calls.append({
        "Pan": n, "Palette": r.Palette,
        "Talc": "Detected" if r["3674 talc OH"] >= 20 and r["667 talc Mg-OH"] > 0 else "Not detected",
        "Mica": ("Detected" if r["mica ratio 3626/3674"] >= 0.40 else
                 "Possible" if r["mica ratio 3626/3674"] >= 0.20 else "Not detected"),
        "Kaolin": "Possible" if r["3689 kaolin OH"] >= 20 else "Not detected",
        "PDMS (dimethicone)": lvl(r["1260 Si-CH3"], 40, 20, 10),
        "Ester": lvl(r["1740 C=O ester"], 20, 10),
        "Metal stearate": lvl(max(r["1540 COO- Zn"], r["1576 COO- Mg"]), 20, 10),
        # Zn 1535-1540 single band; Mg 1573-1578 (SDBS 12704); Ca 1577 + 1540
        "Stearate type": ("" if max(r["1540 COO- Zn"], r["1576 COO- Mg"]) < 10 else
                          "Zn-type" if r["1576 COO- Mg"] < 10 else
                          "Ca-type" if r["1540 COO- Zn"] >= 10 else "Mg-type"),
        "MgCO3": "Detected" if r["747 CO3 nu4/Si-O-Al"] >= 60 and r["1446 CO3 nu3"] >= 10 else "Not detected",
        "CaCO3": "Detected" if r["712 CaCO3 nu4"] >= 25 and r["872 CO3 nu2"] >= 60 else "Not detected",
        "Iron oxide": "Not decidable by FTIR rules (use LIBS)",
    })
C = pd.DataFrame(calls).set_index("Pan")

# ---------------- 3. exploratory PCA / HCA ----------------
def rubber(x, y):
    h = ConvexHull(np.c_[x, y]); v = np.roll(h.vertices, -h.vertices.argmin())
    v = v[: v.argmax() + 1]; v = np.sort(v)
    return np.interp(x, x[v], y[v])

X = []
for n in PANS:
    x, A = load(n)
    a = A - rubber(x, A)
    keep = ~((x > 1880) & (x < 2340))
    a = a[keep]; a = (a - a.mean()) / a.std()          # SNV
    X.append(a)
xk = x[keep]; X = np.array(X); Xc = X - X.mean(0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
ev = S**2 / np.sum(S**2); T = U * S
k = int(np.searchsorted(np.cumsum(ev), 0.95) + 1)
Z = linkage(T[:, :k], method="ward")

fig = plt.figure(figsize=(10, 7.5))
ax1 = fig.add_subplot(2, 2, 1)
col = {"Palette A (1-8)": "#1f4e79", "Haya Queen (9-15)": "#b23a17"}
for i, n in enumerate(PANS):
    pd_ = C.loc[n, "PDMS (dimethicone)"]
    mk = "o" if pd_ in ("Detected", "Minor") else "^"
    ax1.scatter(T[i, 0], T[i, 1], c=col[PALETTE[n]], marker=mk, s=60)
    ax1.annotate(str(n), (T[i, 0], T[i, 1]), textcoords="offset points", xytext=(5, 4), fontsize=8)
ax1.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)"); ax1.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)")
ax1.axhline(0, c="0.8", lw=0.6); ax1.axvline(0, c="0.8", lw=0.6)
ax1.set_title("(a) PCA scores  (o PDMS detected/minor, ^ not)", fontsize=9)
for p, c in col.items(): ax1.scatter([], [], c=c, label=p)
ax1.legend(fontsize=7, loc="best")

ax2 = fig.add_subplot(2, 2, 2)
dendrogram(Z, labels=[str(p) for p in PANS], ax=ax2, color_threshold=0.7 * Z[:, 2].max())
ax2.set_title(f"(b) HCA, Ward, first {k} PCs", fontsize=9); ax2.set_ylabel("Distance")

ax3 = fig.add_subplot(2, 1, 2)
for j, c in ((0, "#1f4e79"), (1, "#b23a17")):
    ax3.plot(xk, Vt[j], c=c, lw=0.8, label=f"PC{j+1} loading")
for b in (3674, 3626, 2917, 1740, 1576, 1540, 1446, 1260, 872, 799, 747, 667, 525):
    ax3.axvline(b, c="0.75", ls=":", lw=0.6)
    ax3.text(b, ax3.get_ylim()[1] if False else 0, "", fontsize=6)
ax3.set_xlim(4000, 400); ax3.set_xlabel("Wavenumber (cm$^{-1}$)"); ax3.set_ylabel("Loading")
ax3.legend(fontsize=7); ax3.set_title("(c) PCA loadings (2340-1880 cm$^{-1}$ removed)", fontsize=9)
fig.tight_layout()
fig.savefig("./eyeshadow_PCA_HCA.png", dpi=300)

# ---------------- 4. export ----------------
with pd.ExcelWriter("./eyeshadow_ftir_screening.xlsx") as w:
    C.to_excel(w, sheet_name="Rule-based calls")
    F.round(2).to_excel(w, sheet_name="Band S-N features")
    pd.DataFrame({"PC": [f"PC{i+1}" for i in range(6)],
                  "Explained variance (%)": (ev[:6] * 100).round(1)}).to_excel(
        w, sheet_name="PCA variance", index=False)
    pd.DataFrame({"Pan": PANS, **{f"PC{i+1}": T[:, i].round(3) for i in range(3)}}).to_excel(
        w, sheet_name="PCA scores", index=False)
print(C.to_string()); print("EV %", (ev[:5] * 100).round(1), "k =", k)
