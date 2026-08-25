#!/usr/bin/env python3
"""
Composite Figure: External Validation of 37-Gene Signature
===========================================================
"""

import numpy as np
# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use run_figure5.py.
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS
# ============================================================================
MSBB_PRED_LOOCV = "ADKP/results/method_3_predictions.csv"
MSBB_PRED_5FOLD = "ADKP/results/method_4_predictions.csv"
MSBB_COEF_LOOCV = "ADKP/results/method_3_coefficients.csv"

GSE_PRED_LOOCV = "alzheimer_meta/new_meta_analysis/results/method_3_predictions.csv"
GSE_PRED_5FOLD = "alzheimer_meta/new_meta_analysis/results/method_4_predictions.csv"
GSE_COEF_LOOCV = "alzheimer_meta/new_meta_analysis/results/method_3_coefficients.csv"

OUTPUT_JPEG = "alzheimer_meta/ad_meta/FIGURES/figure_external_validation.tiff"
OUTPUT_PDF = "alzheimer_meta/ad_meta/FIGURES/figure_external_validation.pdf"

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data...")
msbb_m3 = pd.read_csv(MSBB_PRED_LOOCV)
msbb_m4 = pd.read_csv(MSBB_PRED_5FOLD)
msbb_coef = pd.read_csv(MSBB_COEF_LOOCV)

gse_m3 = pd.read_csv(GSE_PRED_LOOCV)
gse_m4 = pd.read_csv(GSE_PRED_5FOLD)
gse_coef = pd.read_csv(GSE_COEF_LOOCV)

# ============================================================================
# COMPUTE METRICS
# ============================================================================
def get_metrics(df):
    y_true, y_pred, y_proba = df['y_true'], df['y_pred'], df['y_proba']
    return {
        'auc': roc_auc_score(y_true, y_proba),
        'accuracy': accuracy_score(y_true, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0),
        'f1': f1_score(y_true, y_pred, zero_division=0),
    }

msbb_metrics = get_metrics(msbb_m3)
gse_metrics = get_metrics(gse_m3)

def compute_fold_aucs(df, n_folds=5):
    n = len(df)
    fold_size = n // n_folds
    aucs = []
    for i in range(n_folds):
        start = i * fold_size
        end = start + fold_size if i < n_folds - 1 else n
        fold = df.iloc[start:end]
        try:
            aucs.append(roc_auc_score(fold['y_true'], fold['y_proba']))
        except ValueError:
            aucs.append(np.nan)
    return aucs

msbb_fold_aucs = compute_fold_aucs(msbb_m4)
gse_fold_aucs = compute_fold_aucs(gse_m4)

# ============================================================================
# STYLE — large fonts, small canvas so text appears big relative to figure
# ============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 16,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': '#333333',
    'axes.linewidth': 1.2,
})

MSBB_COLOR = '#2C73B4'
GSE_COLOR = '#E85D2A'

# ============================================================================
# FIGURE — compact canvas so fonts appear large
# ============================================================================
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.55, wspace=0.45,
                       left=0.09, right=0.97, top=0.90, bottom=0.07)

# ---- helper: panel letter top-left, legend top-right above axes ----
def add_panel_letter(ax, letter):
    """Place panel letter at top-left, well above the axes."""
    ax.text(-0.12, 1.15, letter, transform=ax.transAxes,
            fontsize=26, fontweight='bold', va='top', ha='left')

def legend_above_right(ax, ncol=1, fontsize=13):
    """Place legend above axes, anchored to top-right so it doesn't clash with panel letter."""
    ax.legend(loc='lower right', bbox_to_anchor=(1.0, 1.02), ncol=ncol,
              frameon=True, fancybox=True, framealpha=0.9, fontsize=fontsize)

# ------ Panel A: ROC Curves ------
ax_a = fig.add_subplot(gs[0, 0])

msbb_fpr, msbb_tpr, _ = roc_curve(msbb_m3['y_true'], msbb_m3['y_proba'])
gse_fpr, gse_tpr, _ = roc_curve(gse_m3['y_true'], gse_m3['y_proba'])

ax_a.plot(msbb_fpr, msbb_tpr, color=MSBB_COLOR, linewidth=2.5,
          label=f'MSBB (AUC = {msbb_metrics["auc"]:.3f})', zorder=3)
ax_a.plot(gse_fpr, gse_tpr, color=GSE_COLOR, linewidth=2.5,
          label=f'GSE125583 (AUC = {gse_metrics["auc"]:.3f})', zorder=3)
ax_a.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.4, label='Random (AUC = 0.5)')
ax_a.fill_between(msbb_fpr, msbb_tpr, alpha=0.08, color=MSBB_COLOR)
ax_a.fill_between(gse_fpr, gse_tpr, alpha=0.08, color=GSE_COLOR)

ax_a.set_xlabel('False Positive Rate')
ax_a.set_ylabel('True Positive Rate')
add_panel_letter(ax_a, '(A)')
legend_above_right(ax_a, ncol=1, fontsize=11)
ax_a.set_xlim([-0.02, 1.02])
ax_a.set_ylim([-0.02, 1.02])
ax_a.grid(alpha=0.2)
ax_a.set_aspect('equal')

# ------ Panel B: Performance Metrics ------
ax_b = fig.add_subplot(gs[0, 1])

metric_keys = ['auc', 'accuracy', 'balanced_accuracy', 'precision', 'recall', 'f1']
metric_labels = ['AUC', 'Acc', 'Bal Acc', 'Prec', 'Recall', 'F1']
msbb_vals = [msbb_metrics[k] for k in metric_keys]
gse_vals = [gse_metrics[k] for k in metric_keys]

x = np.arange(len(metric_labels))
width = 0.32

bars1 = ax_b.bar(x - width/2, msbb_vals, width, color=MSBB_COLOR, alpha=0.85,
                 label=f'MSBB (n={len(msbb_m3):,})', edgecolor='white', linewidth=0.5, zorder=3)
bars2 = ax_b.bar(x + width/2, gse_vals, width, color=GSE_COLOR, alpha=0.85,
                 label=f'GSE125583 (n={len(gse_m3):,})', edgecolor='white', linewidth=0.5, zorder=3)

for bar, val in zip(bars1, msbb_vals):
    ax_b.text(bar.get_x() + bar.get_width()/2, val + 0.01, f'{val:.2f}',
              ha='center', va='bottom', fontsize=10, fontweight='bold', color=MSBB_COLOR)
for bar, val in zip(bars2, gse_vals):
    ax_b.text(bar.get_x() + bar.get_width()/2, val + 0.01, f'{val:.2f}',
              ha='center', va='bottom', fontsize=10, fontweight='bold', color=GSE_COLOR)

ax_b.axhline(y=0.5, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
ax_b.set_xticks(x)
ax_b.set_xticklabels(metric_labels, rotation=25, ha='right')
ax_b.set_ylabel('Score')
add_panel_letter(ax_b, '(B)')
ax_b.set_ylim([0, 1.12])
legend_above_right(ax_b, ncol=1, fontsize=11)
ax_b.grid(axis='y', alpha=0.2)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ------ Panel C: Fold-wise AUC (BOTTOM-LEFT) ------
ax_c = fig.add_subplot(gs[1, 0])

folds = np.arange(1, 6)
bar_width = 0.32

bars_m = ax_c.bar(folds - bar_width/2, msbb_fold_aucs, bar_width, color=MSBB_COLOR,
                  alpha=0.85, label='MSBB', edgecolor='white', linewidth=0.5, zorder=3)
bars_g = ax_c.bar(folds + bar_width/2, gse_fold_aucs, bar_width, color=GSE_COLOR,
                  alpha=0.85, label='GSE125583', edgecolor='white', linewidth=0.5, zorder=3)

ax_c.axhline(y=np.nanmean(msbb_fold_aucs), color=MSBB_COLOR, linestyle='-', linewidth=1.8,
             alpha=0.7, label=f'MSBB mean: {np.nanmean(msbb_fold_aucs):.3f}')
ax_c.axhline(y=np.nanmean(gse_fold_aucs), color=GSE_COLOR, linestyle='-', linewidth=1.8,
             alpha=0.7, label=f'GSE mean: {np.nanmean(gse_fold_aucs):.3f}')
ax_c.axhline(y=0.5, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)

for bar, val in zip(bars_m, msbb_fold_aucs):
    ax_c.text(bar.get_x() + bar.get_width()/2, val + 0.008, f'{val:.2f}',
              ha='center', va='bottom', fontsize=11, color=MSBB_COLOR, fontweight='bold')
for bar, val in zip(bars_g, gse_fold_aucs):
    ax_c.text(bar.get_x() + bar.get_width()/2, val + 0.008, f'{val:.2f}',
              ha='center', va='bottom', fontsize=11, color=GSE_COLOR, fontweight='bold')

ax_c.set_xlabel('Fold')
ax_c.set_ylabel('AUC')
add_panel_letter(ax_c, '(C)')
ax_c.set_xticks(folds)
# FIX: just numbers 1–5 instead of "Fold 1" etc.
ax_c.set_xticklabels(['1', '2', '3', '4', '5'])
ax_c.set_ylim([0.45, 1.05])
legend_above_right(ax_c, ncol=2, fontsize=10)
ax_c.grid(axis='y', alpha=0.2)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ------ Panel D: Coefficient Comparison / Dumbbell Plot (BOTTOM-RIGHT) ------
ax_d = fig.add_subplot(gs[1, 1])

msbb_coef_dict = dict(zip(msbb_coef['gene'], msbb_coef['mean_coefficient']))
gse_coef_dict = dict(zip(gse_coef['gene'], gse_coef['mean_coefficient']))

top_msbb = sorted(msbb_coef_dict.keys(), key=lambda g: abs(msbb_coef_dict[g]), reverse=True)[:12]
top_gse = sorted(gse_coef_dict.keys(), key=lambda g: abs(gse_coef_dict[g]), reverse=True)[:12]
all_top = list(dict.fromkeys(top_msbb + top_gse))[:15]
all_top.sort(key=lambda g: (abs(msbb_coef_dict.get(g, 0)) + abs(gse_coef_dict.get(g, 0))) / 2)

y_pos = np.arange(len(all_top))
msbb_c = [msbb_coef_dict.get(g, 0) for g in all_top]
gse_c = [gse_coef_dict.get(g, 0) for g in all_top]

for i, (m, g) in enumerate(zip(msbb_c, gse_c)):
    ax_d.plot([m, g], [i, i], color='#CCCCCC', linewidth=1.5, zorder=1)

ax_d.scatter(msbb_c, y_pos, color=MSBB_COLOR, s=100, zorder=3,
             edgecolors='white', linewidth=0.7, label='MSBB')
ax_d.scatter(gse_c, y_pos, color=GSE_COLOR, s=100, zorder=3,
             edgecolors='white', linewidth=0.7, marker='D', label='GSE125583')

ax_d.axvline(x=0, color='black', linewidth=0.8, zorder=2)
ax_d.set_yticks(y_pos)
ax_d.set_yticklabels(all_top, fontsize=13)
ax_d.set_xlabel('Elastic-Net Coefficient')
add_panel_letter(ax_d, '(D)')
legend_above_right(ax_d, ncol=2, fontsize=12)
ax_d.grid(axis='x', alpha=0.2)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

xlim = ax_d.get_xlim()
ax_d.axvspan(0, xlim[1], alpha=0.05, color='red')
ax_d.axvspan(xlim[0], 0, alpha=0.05, color='blue')
ax_d.text(xlim[1] * 0.75, 0.5, 'AD ↑', fontsize=14,
          color='red', alpha=0.5, ha='center', fontweight='bold')
ax_d.text(xlim[0] * 0.75, 0.5, 'Control ↑', fontsize=14,
          color='blue', alpha=0.5, ha='center', fontweight='bold')
ax_d.set_xlim(xlim)

# ------ Save ------
plt.savefig(OUTPUT_JPEG, dpi=900, bbox_inches='tight', facecolor='white',
            format='tiff', pil_kwargs={'quality': 95})
plt.savefig(OUTPUT_PDF, dpi=900, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved: {OUTPUT_JPEG}")
print(f"✓ Saved: {OUTPUT_PDF}")
plt.show()
