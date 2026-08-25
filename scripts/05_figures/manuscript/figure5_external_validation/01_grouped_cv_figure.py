#!/usr/bin/env python3
"""
Figure 5: Validation performance of the 37-gene signature in independent cohorts
===============================================================================
All four panels now show DONOR-GROUPED 5-fold cross-validation, so the figure is
internally consistent (minor comment 3). Panel B adds sensitivity and specificity
with bootstrap 95% CIs (minor comment 4).

CHANGES FROM THE PREVIOUS VERSION
  - Panels A/B previously showed LOOCV (method_3) while panel C showed 5-fold
    (method_4). All panels are now the same analysis.
  - compute_fold_aucs() previously sliced the predictions file into 5 contiguous
    chunks and treated them as folds. The predictions file is in SAMPLE order,
    not fold order, so those were not the real folds. Fold AUCs are now taken
    from the CV loop's own output (FOLD_AUCS below) or from a 'fold' column.
  - Precision/Recall replaced by Sensitivity/Specificity in panel B.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix
import warnings
from pathlib import Path
import json
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS  -- point these at the grouped_cv outputs
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RES = PROJECT_ROOT / "results" / "validation" / "grouped_cv"

MSBB_PRED  = RES / "MSBB_grouped_cv_predictions.csv"
MSBB_COEF  = RES / "MSBB_grouped_cv_coefficients.csv"
MSBB_MET   = RES / "MSBB_grouped_cv_metrics.csv"

GSE_PRED   = RES / "GSE125583_grouped_cv_predictions.csv"
GSE_COEF   = RES / "GSE125583_grouped_cv_coefficients.csv"
GSE_MET    = RES / "GSE125583_grouped_cv_metrics.csv"

OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "main"
OUTPUT_TIFF = OUTPUT_DIR / "Figure5_external_validation.tiff"
OUTPUT_PDF  = OUTPUT_DIR / "Figure5_external_validation.pdf"
OUTPUT_PNG  = OUTPUT_DIR / "Figure5_external_validation.png"

def read_manifest(cohort):
    with open(RES / f"{cohort}_grouped_cv_manifest.json") as handle:
        return json.load(handle)

# ============================================================================
# LOAD
# ============================================================================
print("Loading grouped-CV outputs...")
msbb_pred = pd.read_csv(MSBB_PRED)
gse_pred  = pd.read_csv(GSE_PRED)
msbb_manifest = read_manifest("MSBB")
gse_manifest = read_manifest("GSE125583")


def read_coef(path):
    """Coefficient files are written from a Series, so the gene column may be
    unnamed. Normalise to columns ['gene', 'coef']."""
    d = pd.read_csv(path)
    d.columns = ['gene', 'coef'] + list(d.columns[2:])
    return dict(zip(d['gene'], d['coef']))


msbb_coef_dict = read_coef(MSBB_COEF)
gse_coef_dict  = read_coef(GSE_COEF)


def read_ci(path):
    """metrics.csv -> {metric: (value, lo, hi)} keyed by lowercase name."""
    try:
        d = pd.read_csv(path)
    except FileNotFoundError:
        return {}
    return {r['metric'].strip().lower(): (r['value'], r['ci_low'], r['ci_high'])
            for _, r in d.iterrows()}


msbb_ci = read_ci(MSBB_MET)
gse_ci  = read_ci(GSE_MET)


def get_metrics(df):
    y, p, pr = df['y_true'].values, df['y_pred'].values, df['y_proba'].values
    tn, fp, fn, tp = confusion_matrix(y, p, labels=[0, 1]).ravel()
    return {
        'auc':               roc_auc_score(y, pr),
        'accuracy':          (tp + tn) / len(y),
        'balanced accuracy': 0.5 * (tp / (tp + fn) + tn / (tn + fp)),
        'sensitivity':       tp / (tp + fn),
        'specificity':       tn / (tn + fp),
        'f1':                2 * tp / (2 * tp + fp + fn),
    }


msbb_metrics = get_metrics(msbb_pred)
gse_metrics  = get_metrics(gse_pred)


def fold_aucs(cohort, df):
    """Compute fold AUCs from actual CV assignments written by the model."""
    if 'fold' not in df.columns:
        raise SystemExit(f"{cohort} predictions lack the required 'fold' column")
    return [roc_auc_score(g['y_true'], g['y_proba'])
            for _, g in df.groupby('fold', sort=True)]


msbb_fold = fold_aucs('MSBB', msbb_pred)
gse_fold  = fold_aucs('GSE125583', gse_pred)

# ============================================================================
# STYLE
# ============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10, 'axes.titlesize': 10, 'axes.titleweight': 'bold',
    'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 8, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'axes.edgecolor': '#333333', 'axes.linewidth': 1.0,
})
MSBB_COLOR, GSE_COLOR = '#2C73B4', '#E85D2A'

fig = plt.figure(figsize=(7, 7))
gs = gridspec.GridSpec(2, 2, hspace=0.55, wspace=0.45,
                       left=0.09, right=0.97, top=0.90, bottom=0.07)


def panel_letter(ax, s):
    ax.text(-0.12, 1.15, s, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')


def legend_above(ax, ncol=1, fontsize=8):
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), ncol=ncol,
              frameon=True, fancybox=True, framealpha=0.9, fontsize=fontsize)


# ------ Panel A: ROC ------
ax = fig.add_subplot(gs[0, 0])
for df, col, lab, m in [(msbb_pred, MSBB_COLOR, 'MSBB', msbb_metrics),
                        (gse_pred, GSE_COLOR, 'GSE125583', gse_metrics)]:
    fpr, tpr, _ = roc_curve(df['y_true'], df['y_proba'])
    ax.plot(fpr, tpr, color=col, linewidth=1.8,
            label=f'{lab} (AUC = {m["auc"]:.3f})', zorder=3)
    ax.fill_between(fpr, tpr, alpha=0.08, color=col)
ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.4, label='Random (AUC = 0.5)')
ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
panel_letter(ax, '(A)'); legend_above(ax, fontsize=7)
ax.set_xlim([-0.02, 1.02]); ax.set_ylim([-0.02, 1.02])
ax.grid(alpha=0.2); ax.set_aspect('equal')

# ------ Panel B: metrics with bootstrap 95% CI ------
ax = fig.add_subplot(gs[0, 1])
keys   = ['auc', 'accuracy', 'balanced accuracy', 'sensitivity', 'specificity', 'f1']
labels = ['AUC', 'Acc', 'Bal Acc', 'Sens', 'Spec', 'F1']
x, width = np.arange(len(labels)), 0.32


def err_bars(ci, met):
    """Asymmetric error bars from the metrics CSV; zeros where CI unavailable."""
    lo, hi = [], []
    for k in keys:
        if k in ci:
            v, l, h = ci[k]
            lo.append(max(v - l, 0)); hi.append(max(h - v, 0))
        else:
            lo.append(0); hi.append(0)
    return np.array([lo, hi])


for vals, ci, off, col, lab, n in [
        ([msbb_metrics[k] for k in keys], msbb_ci, -width/2, MSBB_COLOR,
         f'MSBB ({msbb_manifest["n_samples"]:,} samples, {msbb_manifest["n_donors"]} donors)', None),
        ([gse_metrics[k] for k in keys], gse_ci, width/2, GSE_COLOR,
         f'GSE125583 (n={len(gse_pred):,})', None)]:
    b = ax.bar(x + off, vals, width, color=col, alpha=0.85, label=lab,
               edgecolor='white', linewidth=0.5, zorder=3,
               yerr=err_bars(ci, vals), capsize=2,
               error_kw={'elinewidth': 0.8, 'ecolor': '#444444'})
    for bar, val in zip(b, vals):
        hi = err_bars(ci, vals)[1][list(b).index(bar)]
        ax.text(bar.get_x() + bar.get_width()/2, val + hi + 0.025, f'{val:.2f}',
                ha='center', va='bottom', fontsize=6.5, fontweight='bold', color=col)

ax.axhline(0.5, color='grey', linestyle='--', linewidth=0.7, alpha=0.5)
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=25, ha='right')
ax.set_ylabel('Score'); panel_letter(ax, '(B)')
ax.set_ylim([0, 1.15]); legend_above(ax, fontsize=6.5)
ax.grid(axis='y', alpha=0.2)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ------ Panel C: fold-wise AUC ------
ax = fig.add_subplot(gs[1, 0])
folds, bw = np.arange(1, 6), 0.32
bm = ax.bar(folds - bw/2, msbb_fold, bw, color=MSBB_COLOR, alpha=0.85,
            label='MSBB', edgecolor='white', linewidth=0.5, zorder=3)
bg = ax.bar(folds + bw/2, gse_fold, bw, color=GSE_COLOR, alpha=0.85,
            label='GSE125583', edgecolor='white', linewidth=0.5, zorder=3)
ax.axhline(np.mean(msbb_fold), color=MSBB_COLOR, linewidth=1.3, alpha=0.7,
           label=f'MSBB mean: {np.mean(msbb_fold):.3f}')
ax.axhline(np.mean(gse_fold), color=GSE_COLOR, linewidth=1.3, alpha=0.7,
           label=f'GSE mean: {np.mean(gse_fold):.3f}')
ax.axhline(0.5, color='grey', linestyle='--', linewidth=0.7, alpha=0.5)
for bars, vals, col in [(bm, msbb_fold, MSBB_COLOR), (bg, gse_fold, GSE_COLOR)]:
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.008, f'{val:.2f}',
                ha='center', va='bottom', fontsize=7, color=col, fontweight='bold')
ax.set_xlabel('Fold'); ax.set_ylabel('AUC'); panel_letter(ax, '(C)')
ax.set_xticks(folds); ax.set_ylim([0.45, 1.05])
ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), ncol=2, frameon=True,
          fancybox=True, framealpha=0.9, fontsize=7,
          columnspacing=1.0, handletextpad=0.4)
ax.grid(axis='y', alpha=0.2)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ------ Panel D: coefficient dumbbell ------
ax = fig.add_subplot(gs[1, 1])
top_m = sorted(msbb_coef_dict, key=lambda g: abs(msbb_coef_dict[g]), reverse=True)[:12]
top_g = sorted(gse_coef_dict,  key=lambda g: abs(gse_coef_dict[g]),  reverse=True)[:12]
genes = list(dict.fromkeys(top_m + top_g))[:15]
genes.sort(key=lambda g: (abs(msbb_coef_dict.get(g, 0)) + abs(gse_coef_dict.get(g, 0))) / 2)

y_pos = np.arange(len(genes))
mc = [msbb_coef_dict.get(g, 0) for g in genes]
gc = [gse_coef_dict.get(g, 0) for g in genes]
for i, (m, g) in enumerate(zip(mc, gc)):
    ax.plot([m, g], [i, i], color='#CCCCCC', linewidth=1.2, zorder=1)
ax.scatter(mc, y_pos, color=MSBB_COLOR, s=55, zorder=3,
           edgecolors='white', linewidth=0.7, label='MSBB')
ax.scatter(gc, y_pos, color=GSE_COLOR, s=55, zorder=3,
           edgecolors='white', linewidth=0.7, marker='D', label='GSE125583')
ax.axvline(0, color='black', linewidth=0.7, zorder=2)
ax.set_yticks(y_pos); ax.set_yticklabels(genes, fontsize=8)
ax.set_xlabel('Elastic-Net Coefficient'); panel_letter(ax, '(D)')
legend_above(ax, ncol=2, fontsize=7)
ax.grid(axis='x', alpha=0.2)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
xl = ax.get_xlim()
ax.axvspan(0, xl[1], alpha=0.05, color='red')
ax.axvspan(xl[0], 0, alpha=0.05, color='blue')
ax.text(xl[1]*0.75, 0.5, 'AD ↑', fontsize=9, color='red', alpha=0.5,
        ha='center', fontweight='bold')
ax.text(xl[0]*0.75, 0.5, 'Control ↑', fontsize=9, color='blue', alpha=0.5,
        ha='center', fontweight='bold')
ax.set_xlim(xl)

# ------ Save ------
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_TIFF, dpi=900, bbox_inches='tight', facecolor='white',
            format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
plt.savefig(OUTPUT_PDF, dpi=900, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved {OUTPUT_TIFF}\nSaved {OUTPUT_PDF}\nSaved {OUTPUT_PNG}")
