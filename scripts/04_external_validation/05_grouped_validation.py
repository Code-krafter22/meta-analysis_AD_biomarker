#!/usr/bin/env python3
"""
Donor-grouped validation of the 37-gene signature in MSBB
=========================================================

Addresses Reviewer 2, major comment 3:
  - MSBB samples are grouped by donor (individualID), so all samples from a
    donor fall in the same CV fold. Previously StratifiedKFold split at the
    SAMPLE level, placing e.g. a donor's BM22 sample in training and their
    BM36 sample in test -> information leakage -> optimistic AUC.
  - Confidence intervals are computed by resampling DONORS, not samples.
    Resampling samples treats 4 correlated tissues from one donor as 4
    independent observations and gives intervals that are too narrow.
  - MODE 'external' trains the full model in one cohort and applies it to the
    other with no retraining. That is genuine external validation.

MODES
  'grouped_cv'  : donor-grouped stratified 5-fold CV within a cohort  [primary]
  'logo'        : leave-one-donor-out CV                              [sensitivity]
  'external'    : train on cohort A, apply unchanged to cohort B      [primary external]

Hyperparameters are FIXED a priori (L1_RATIO, C_VALUE). No grid search anywhere,
which is what Methods 2.6 claims. Do not reintroduce LogisticRegressionCV.
"""

import numpy as np
import pandas as pd
import os
import json
import platform
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_auc_score, accuracy_score, balanced_accuracy_score,
    precision_score, recall_score, f1_score, confusion_matrix, roc_curve
)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
MODE = os.getenv('VALIDATION_MODE', 'grouped_cv')
COHORT = os.getenv('VALIDATION_COHORT', 'GSE125583')
EXTERNAL_DIRECTION = os.getenv('EXTERNAL_DIRECTION', 'GSE125583->MSBB')

L1_RATIO = 0.5
C_VALUE = 1.0
N_FOLDS = 5
N_BOOT = 1000
RANDOM_STATE = 42

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MSBB_EXPR = os.getenv('MSBB_EXPR')
MSBB_META = os.getenv('MSBB_META')
GSE_EXPR = PROJECT_ROOT / 'data/processed/GSE125583/GSE125583_log2cpm_SELECTED_GENES.csv'
GSE_META = PROJECT_ROOT / 'data/metadata/GSE125583/metadata_200samples.csv'
OUTDIR = Path(os.getenv('VALIDATION_OUTPUT_DIR',
                         PROJECT_ROOT / 'results/validation/grouped_cv'))

np.random.seed(RANDOM_STATE)

GENES_37 = [
    "ADAM33","AEBP1","CCDC102A","CHML","CLDN9","COL25A1","ELOVL4","GAD1","GAD2",
    "GFAP","HPRT1","HSPB1","HSPB7","ITFG1","KANK2","KLF15","MAS1","MRGPRF",
    "NAP1L5","NCALD","NEUROD6","NRN1","NUPR1","OPN3","PIK3R5","PRELP","PRX",
    "RAB3B","RAB3C","RGS4","RPH3A","SCG2","SERPINI1","STAT4","TCEA3","TMPRSS5",
    "TRIM36"
]
# This is the definitive 37-gene vector used for validation. GSE125583 contains
# 36/37 genes because COL25A1 is unavailable; MSBB contains all 37.


def load_expression(path):
    """Load expression, orient to samples x genes."""
    d = pd.read_csv(path, index_col=0)
    # genes-as-rows if few rows and many columns
    if d.shape[0] < 100 and d.shape[1] > 100:
        d = d.T
    return d


def load_msbb():
    """MSBB: returns X (samples x genes), y (0/1), groups (donor id)."""
    if not MSBB_EXPR or not MSBB_META:
        raise SystemExit(
            "MSBB is controlled-access. Set MSBB_EXPR and MSBB_META to local "
            "authorized files; see scripts/04_external_validation/README.md."
        )
    expr = load_expression(MSBB_EXPR)
    meta = pd.read_csv(MSBB_META)

    # --- CHANGE 3: filter BEFORE any splitting -----------------------------
    # 'Uncertain' donors (n=227 samples) and flagged exclusions must be removed
    # here, not after the split, or they leak into fold construction.
    meta = meta[meta['Diagnosis'].isin(['AD', 'Control'])]
    if 'exclude' in meta.columns:
        excl = meta['exclude'].astype(str).str.upper().isin(['TRUE', '1'])
        meta = meta[~excl]

    meta = meta.set_index('specimenID')
    common = expr.index.intersection(meta.index)
    if len(common) == 0:
        raise SystemExit(
            "No overlap between expression columns and specimenID.\n"
            f"  expression ids: {list(expr.index[:3])}\n"
            f"  metadata ids:   {list(meta.index[:3])}\n"
            "MSBB ids sometimes carry suffixes - strip them before joining."
        )
    expr, meta = expr.loc[common], meta.loc[common]

    genes = [g for g in GENES_37 if g in expr.columns]
    X = expr[genes].copy()
    y = (meta['Diagnosis'] == 'AD').astype(int)
    groups = meta['individualID'].values      # <-- the grouping variable

    ok = ~X.isna().any(axis=1)
    X, y, groups = X[ok], y[ok], groups[ok.values]

    print(f"  MSBB: {X.shape[0]} samples from {pd.Series(groups).nunique()} donors, "
          f"{X.shape[1]}/{len(GENES_37)} genes")
    print(f"        AD {int((y==1).sum())}, Control {int((y==0).sum())}")
    spd = pd.Series(groups).value_counts()
    print(f"        samples/donor: median {int(spd.median())}, max {int(spd.max())}")
    return X, y.values, groups


def load_gse():
    """GSE125583: one sample per donor, so groups == sample index."""
    expr = load_expression(GSE_EXPR)
    meta = pd.read_csv(GSE_META, index_col=0)
    common = expr.index.intersection(meta.index)
    expr, meta = expr.loc[common], meta.loc[common]

    dcol = 'diagnosis:ch1' if 'diagnosis:ch1' in meta.columns else meta.columns[0]
    lab = meta[dcol].astype(str).str.lower()
    y = pd.Series(np.where(lab.str.contains('alzh|^ad$|disease'), 1,
                  np.where(lab.str.contains('control|normal|^cn$'), 0, -1)),
                  index=meta.index)

    genes = [g for g in GENES_37 if g in expr.columns]
    X = expr[genes].copy()
    keep = (y >= 0) & ~X.isna().any(axis=1)
    X, y = X[keep], y[keep]
    groups = np.asarray(X.index)   # one sample per donor

    print(f"  GSE125583: {X.shape[0]} samples, {X.shape[1]}/{len(GENES_37)} genes")
    print(f"             AD {int((y==1).sum())}, Control {int((y==0).sum())}")
    return X, y.values, groups


def make_pipeline():
    """Scaling nested inside the pipeline so it is refit per training fold."""
    return Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegression(
            penalty='elasticnet', solver='saga',
            l1_ratio=L1_RATIO, C=C_VALUE,          # fixed a priori - no tuning
            class_weight='balanced',
            max_iter=10000, random_state=RANDOM_STATE))
    ])


def metrics(y_true, y_pred, y_proba):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        'AUC': roc_auc_score(y_true, y_proba),
        'Accuracy': accuracy_score(y_true, y_pred),
        'Balanced accuracy': balanced_accuracy_score(y_true, y_pred),
        'Sensitivity': tp / (tp + fn) if (tp + fn) else np.nan,
        'Specificity': tn / (tn + fp) if (tn + fp) else np.nan,
        'Precision': precision_score(y_true, y_pred, zero_division=0),
        'F1': f1_score(y_true, y_pred, zero_division=0),
    }


def donor_bootstrap(y_true, y_pred, y_proba, groups, n_boot=N_BOOT, seed=0):
    """
    CHANGE 4: resample DONORS with replacement, then take every sample
    belonging to each drawn donor. Sample-level resampling would treat one
    donor's four regions as four independent observations and understate the
    interval width.
    """
    rng = np.random.default_rng(seed)
    donors = np.unique(groups)
    idx_by_donor = {d: np.where(groups == d)[0] for d in donors}
    keys = list(metrics(y_true, y_pred, y_proba).keys())
    boot = {k: [] for k in keys}

    for _ in range(n_boot):
        drawn = rng.choice(donors, size=len(donors), replace=True)
        idx = np.concatenate([idx_by_donor[d] for d in drawn])
        if len(np.unique(y_true[idx])) < 2:
            continue
        m = metrics(y_true[idx], y_pred[idx], y_proba[idx])
        for k in keys:
            boot[k].append(m[k])

    return {k: (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5))
            for k, v in boot.items() if len(v)}


def report(name, y_true, y_pred, y_proba, groups):
    m = metrics(y_true, y_pred, y_proba)
    ci = donor_bootstrap(y_true, y_pred, y_proba, groups)
    print(f"\n{name}")
    print("-" * 62)
    rows = []
    for k, v in m.items():
        lo, hi = ci.get(k, (np.nan, np.nan))
        print(f"  {k:<20s} {v:.3f}   95% CI {lo:.3f}-{hi:.3f}")
        rows.append({'metric': k, 'value': v, 'ci_low': lo, 'ci_high': hi})
    return pd.DataFrame(rows)


# ============================================================================
# MODES
# ============================================================================
def run_grouped_cv(X, y, groups, label):
    """CHANGE 1: StratifiedGroupKFold, and pass groups to .split()."""
    print(f"\nDonor-grouped stratified {N_FOLDS}-fold CV  [{label}]")
    cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True,
                              random_state=RANDOM_STATE)

    proba = np.zeros(len(y)); pred = np.zeros(len(y), dtype=int)
    fold_id = np.zeros(len(y), dtype=int)
    coefs, fold_auc = [], []

    # NOTE the third argument. Omitting `groups` silently reverts to
    # ungrouped behaviour and reintroduces the leak.
    for i, (tr, te) in enumerate(cv.split(X, y, groups), 1):
        assert not (set(groups[tr]) & set(groups[te])), "donor in both folds"
        pipe = make_pipeline().fit(X.iloc[tr], y[tr])
        proba[te] = pipe.predict_proba(X.iloc[te])[:, 1]
        pred[te] = pipe.predict(X.iloc[te])
        fold_id[te] = i
        a = roc_auc_score(y[te], proba[te]); fold_auc.append(a)
        print(f"  fold {i}: AUC {a:.3f}  "
              f"({len(tr)} train / {len(te)} test samples, "
              f"{len(set(groups[te]))} held-out donors)")
        coefs.append(pd.Series(pipe.named_steps['clf'].coef_[0], index=X.columns))

    print(f"  mean fold AUC {np.mean(fold_auc):.3f} +/- {np.std(fold_auc):.3f}")
    df = report(f"Pooled out-of-fold performance [{label}]", y, pred, proba, groups)
    coef_df = pd.DataFrame(coefs).mean().sort_values(key=abs, ascending=False)
    return df, coef_df, proba, pred, fold_id, fold_auc


def run_logo(X, y, groups, label):
    """Leave-one-DONOR-out. Leave-one-sample-out leaks the same donor."""
    print(f"\nLeave-one-donor-out CV  [{label}]  "
          f"({pd.Series(groups).nunique()} iterations)")
    logo = LeaveOneGroupOut()
    proba = np.zeros(len(y)); pred = np.zeros(len(y), dtype=int)
    for j, (tr, te) in enumerate(logo.split(X, y, groups), 1):
        pipe = make_pipeline().fit(X.iloc[tr], y[tr])
        proba[te] = pipe.predict_proba(X.iloc[te])[:, 1]
        pred[te] = pipe.predict(X.iloc[te])
        if j % 50 == 0:
            print(f"  {j} donors done", end='\r')
    print()
    return report(f"Leave-one-donor-out [{label}]", y, pred, proba, groups)


def run_external(direction):
    """
    CHANGE 5 / M3 part 3: train the complete model in one cohort, apply it to
    the other with NO retraining. The scaler is fit on the training cohort only
    and applied unchanged - re-standardising on the test cohort would leak.
    """
    print(f"\nFully external validation: {direction}")
    Xm, ym, gm = load_msbb()
    Xg, yg, gg = load_gse()

    common = [g for g in Xm.columns if g in Xg.columns]
    print(f"  genes common to both cohorts: {len(common)}/{len(GENES_37)}")
    missing = sorted(set(GENES_37) - set(common))
    if missing:
        print(f"  not available in both: {', '.join(missing)}")

    if direction == 'MSBB->GSE125583':
        Xtr, ytr, Xte, yte, gte, lab = Xm[common], ym, Xg[common], yg, gg, 'GSE125583'
    else:
        Xtr, ytr, Xte, yte, gte, lab = Xg[common], yg, Xm[common], ym, gm, 'MSBB'

    pipe = make_pipeline().fit(Xtr, ytr)          # trained once, never refit
    proba = pipe.predict_proba(Xte)[:, 1]
    pred = pipe.predict(Xte)

    df = report(f"External: trained on {direction.split('->')[0]}, "
                f"applied to {lab}", yte, pred, proba, gte)
    coef = pd.Series(pipe.named_steps['clf'].coef_[0],
                     index=Xtr.columns).sort_values(key=abs, ascending=False)
    return df, coef, proba, pred, yte


# ============================================================================
# MAIN
# ============================================================================
if __name__ == '__main__':
    os.makedirs(OUTDIR, exist_ok=True)
    print("=" * 62)
    print(f"MODE: {MODE}")
    print("=" * 62)

    if MODE == 'external':
        df, coef, proba, pred, ytrue = run_external(EXTERNAL_DIRECTION)
        tag = 'external_' + EXTERNAL_DIRECTION.replace('->', '_to_')
        coef.to_frame('coefficient').to_csv(f'{OUTDIR}/{tag}_coefficients.csv')
        pd.DataFrame({'y_true': ytrue, 'y_pred': pred,
                      'y_proba': proba}).to_csv(f'{OUTDIR}/{tag}_predictions.csv',
                                                index=False)
    else:
        print("\nLoading data...")
        X, y, groups = load_msbb() if COHORT == 'MSBB' else load_gse()
        tag = f'{COHORT}_{MODE}'
        if MODE == 'grouped_cv':
            df, coef, proba, pred, fold_id, fold_auc = run_grouped_cv(X, y, groups, COHORT)
            coef.to_frame('mean_coefficient').to_csv(f'{OUTDIR}/{tag}_coefficients.csv')
            pd.DataFrame({'y_true': y, 'y_pred': pred,
                          'y_proba': proba, 'fold': fold_id}).to_csv(f'{OUTDIR}/{tag}_predictions.csv',
                                                    index=False)
            with open(f'{OUTDIR}/{tag}_manifest.json', 'w') as handle:
                json.dump({
                    'cohort': COHORT,
                    'n_samples': int(len(y)),
                    'n_donors': int(pd.Series(groups).nunique()),
                    'n_genes': int(X.shape[1]),
                    'n_folds': N_FOLDS,
                    'fold_aucs': [float(v) for v in fold_auc],
                    'l1_ratio': L1_RATIO,
                    'C': C_VALUE,
                    'random_state': RANDOM_STATE,
                    'software': {
                        'python': platform.python_version(),
                        'numpy': np.__version__,
                        'pandas': pd.__version__,
                        'scikit_learn': __import__('sklearn').__version__,
                        'matplotlib': matplotlib.__version__,
                    },
                }, handle, indent=2)
            fpr, tpr, _ = roc_curve(y, proba)
            plt.figure(figsize=(5, 5))
            plt.plot(fpr, tpr, lw=2,
                     label=f'{COHORT} donor-grouped (AUC={roc_auc_score(y, proba):.3f})')
            plt.plot([0, 1], [0, 1], 'k--', lw=1)
            plt.xlabel('False positive rate'); plt.ylabel('True positive rate')
            plt.legend(loc='lower right'); plt.tight_layout()
            plt.savefig(f'{OUTDIR}/{tag}_roc.png', dpi=300)
        else:
            df = run_logo(X, y, groups, COHORT)

    df.to_csv(f'{OUTDIR}/{tag}_metrics.csv', index=False)
    print(f"\nSaved to {OUTDIR}/{tag}_*")
    if COHORT == 'MSBB' and MODE == 'grouped_cv':
        print("\nThe donor-grouped MSBB AUC is expected to be below the prior "
              "sample-level estimate because donor grouping removes leakage.")
