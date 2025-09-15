#!/usr/bin/env python3
import argparse
import json
import os
import joblib
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
# Try to import XGBoost; fall back to scikit-learn if not available (e.g., macOS without libomp)
HAVE_XGB = True
try:
    from xgboost import XGBClassifier
except Exception as e:
    print("XGBoost unavailable, falling back to GradientBoostingClassifier:", e)
    HAVE_XGB = False
    from sklearn.ensemble import GradientBoostingClassifier


NUMERIC_FEATS = [
    'gnomADg_AF', 'heart_lv_tpm', 'heart_aa_tpm', 'GERP++_RS', 'phyloP100way_vertebrate'
]
CAT_FEATS = ['SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred']


def prep(df: pd.DataFrame):
    df = df.copy()
    # Map categorical predictions to numeric
    mapping = {'D':1, 'P':1, 'A':1, 'T':0, 'B':0, 'N':0}
    for c in CAT_FEATS:
        if c in df.columns:
            df[c] = df[c].map(mapping).fillna(0)
    feats = [c for c in NUMERIC_FEATS + CAT_FEATS if c in df.columns]
    X = df[feats].replace([np.inf, -np.inf], np.nan).fillna(0.0).values
    y = df['label'].values.astype(int)
    return X, y, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', required=False)
    ap.add_argument('--data', required=False, help='Alias for --dataset')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    dataset_path = args.dataset or args.data
    if not dataset_path:
        ap.error('one of --dataset or --data is required')

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(dataset_path)
    df = df.dropna(subset=['label'])

    X, y, feats = prep(df)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)

    # Choose model depending on availability
    if HAVE_XGB:
        # A small, fast XGB config
        model = XGBClassifier(
            n_estimators=300,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            n_jobs=4,
            eval_metric='logloss'
        )
        model_type = 'xgboost'
    else:
        # Reasonable sklearn fallback
        model = GradientBoostingClassifier(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.05,
            subsample=0.8 if hasattr(GradientBoostingClassifier(), 'subsample') else 1.0
        )
        model_type = 'sklearn_gradient_boosting'
    model.fit(X_train, y_train)
    proba = model.predict_proba(X_test)[:,1]
    pred = (proba >= 0.5).astype(int)

    metrics = {
        'ROC_AUC': float(roc_auc_score(y_test, proba)),
        'PR_AUC': float(average_precision_score(y_test, proba)),
        'precision': float(precision_score(y_test, pred)),
        'recall': float(recall_score(y_test, pred)),
        'n_test': int(len(y_test)),
        'features': feats,
        'model_type': model_type,
    }
    with open(os.path.join(args.outdir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)

    # Save model
    # Keep filename stable for downstream app regardless of model type
    model_path = os.path.join(args.outdir, 'model_xgb.pkl')
    joblib.dump({'model': model, 'features': feats}, model_path)

    # SHAP
    try:
        import shap
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_test)
        shap_df = pd.DataFrame(shap_values, columns=feats)
        shap_df.head(10).to_csv(os.path.join(args.outdir, 'shap_values.csv'), index=False)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        shap.summary_plot(shap_values, features=X_test, feature_names=feats, show=False)
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, 'shap_summary.png'), dpi=150)
    except Exception as e:
        print('SHAP failed:', e)

    print('Saved model and metrics to', args.outdir)


if __name__ == '__main__':
    main()
