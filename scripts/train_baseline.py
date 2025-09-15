#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features', required=True, help='results/variant_features.csv')
    ap.add_argument('--outdir', required=True, help='results/model')
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.features)
    df = df.dropna(subset=['label'])
    # Minimal feature set
    feats = [c for c in ['gnomADg_AF', 'heart_lv_tpm', 'heart_aa_tpm', 'GERP++_RS', 'phyloP100way_vertebrate'] if c in df.columns]
    # Add simple encodings for SIFT/PolyPhen preds
    for col in ['SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred']:
        if col in df.columns:
            df[col] = df[col].map({'D':1, 'P':1, 'A':1, 'T':0, 'B':0, 'N':0}).fillna(0)
            feats.append(col)

    X = df[feats].fillna(0.0).values
    y = df['label'].values

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    clf = LogisticRegression(max_iter=1000)
    clf.fit(X_train_s, y_train)
    proba = clf.predict_proba(X_test_s)[:,1]
    metrics = {
        'ROC_AUC': float(roc_auc_score(y_test, proba)),
        'PR_AUC': float(average_precision_score(y_test, proba)),
        'n_test': int(len(y_test))
    }
    with open(os.path.join(args.outdir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)
    print('Metrics:', metrics)

    try:
        import shap
        explainer = shap.LinearExplainer(clf, X_train_s, feature_names=feats)
        shap_values = explainer.shap_values(X_test_s)
        shap_df = pd.DataFrame(shap_values, columns=feats)
        shap_out = os.path.join(args.outdir, 'shap_sample.csv')
        shap_df.head(100).to_csv(shap_out, index=False)
        print(f'Wrote SHAP sample to {shap_out}')
    except Exception as e:
        print('SHAP not available or failed:', e)


if __name__ == '__main__':
    main()
