"""260424 Circulant 블로그 포스트용 그림 생성.

두 개의 figure 를 blog/attachments/ 에 저장:
  260424_circulant_01_spectrum.png   — 세 W 의 고유값 크기 |λ_k| 비교
  260424_circulant_02_diffusion.png  — diffusion distance PCA 3D 임베딩
                                        (3 W) × (t=1, 2, 4) grid, color=ring pos

링 60점, 논문 Gaussian / sym cycle / directed shift 세 circulant.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import circulant
from scipy.linalg import orthogonal_procrustes
from sklearn.decomposition import PCA

sys.path.insert(0, "/root/Dropbox/01-rsch/2026-OHS-HST")
from pyhst import l2distance

n = 60
pi = np.pi
ang = np.linspace(-pi, pi - 2 * pi / n, n)
col = plt.get_cmap("rainbow")((ang + pi) / 2 / pi)

vx = np.cos(ang); vy = np.sin(ang)
Σ = l2distance(np.array([vx, vy]).T)
W_gauss = np.exp(-Σ ** 2 / (2 * 0.35 ** 2)) - np.eye(n)

c_cyc = np.zeros(n); c_cyc[1] = c_cyc[n - 1] = 1.0
W_cyc = circulant(c_cyc)

c_shift = np.zeros(n); c_shift[n - 1] = 1.0
W_shift = circulant(c_shift)

Ws = [("Gaussian (paper)", W_gauss),
      ("sym cycle",        W_cyc),
      ("shift (directed)", W_shift)]

out_dir = "/root/Dropbox/01-rsch/2026-OHS-HST/blog/attachments"
os.makedirs(out_dir, exist_ok=True)

# ---------- Figure 1: 스펙트럼 |λ_k| ----------
fig, axes = plt.subplots(1, 3, figsize=(12, 3.2), dpi=150)
for ax, (name, W) in zip(axes, Ws):
    mags = np.abs(np.fft.fft(W[:, 0]))
    ax.stem(np.arange(n), mags, markerfmt="o", basefmt=" ", linefmt="C0-")
    ax.set_xlabel("frequency index $k$")
    ax.set_ylabel(r"$|\hat{c}_k|$")
    ax.set_title(name, fontsize=11)
    ax.grid(True, alpha=0.3)
fig.tight_layout()
path1 = os.path.join(out_dir, "260424_circulant_01_spectrum.png")
fig.savefig(path1, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {path1}")

# ---------- Figure 2: diffusion distance PCA 3D, (W × t) grid ----------
def _pca_align(D, ref=None):
    p = PCA(n_components=3).fit_transform(D)
    if ref is not None:
        R, _ = orthogonal_procrustes(p, ref)
        p = p @ R
    return p

ts = [1, 2, 4]
fig = plt.figure(figsize=(11, 11), dpi=150)
for row, (name, W) in enumerate(Ws):
    rowsum = W.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1
    P = W / rowsum
    ref = None
    for col_idx, t in enumerate(ts):
        Pt = np.linalg.matrix_power(P, t)
        D2 = l2distance(Pt)
        D = np.sqrt(D2)
        p3 = _pca_align(D, ref); ref = p3

        ax = fig.add_subplot(3, 3, row * 3 + col_idx + 1, projection="3d")
        ax.scatter3D(p3[:, 0], p3[:, 1], p3[:, 2], s=22, c=col)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.grid(False)
        uniq = len(np.unique(D2[~np.eye(n, dtype=bool)].round(8)))
        if row == 0:
            ax.set_title(f"$t = {t}$", fontsize=13, pad=8)
        if col_idx == 0:
            ax.text2D(-0.15, 0.5, name, transform=ax.transAxes,
                      ha="right", va="center", fontsize=12,
                      fontstyle="italic", fontweight="bold")
        ax.text2D(0.02, 0.98, f"unique $D^2$: {uniq}",
                  transform=ax.transAxes, ha="left", va="top",
                  fontsize=9, color="#333",
                  bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#aaa", alpha=0.9))

fig.suptitle(r"Diffusion distance PCA 3D embedding — $\mathbf{W}$ × $t$",
             fontsize=14, y=0.995)
fig.tight_layout(rect=(0.03, 0, 1, 0.97))
path2 = os.path.join(out_dir, "260424_circulant_02_diffusion.png")
fig.savefig(path2, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {path2}")
