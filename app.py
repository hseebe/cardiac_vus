import streamlit as st
import pandas as pd
import joblib
import json
import os

st.set_page_config(page_title="Cardiac VUS Demo", layout="wide")

@st.cache_resource
def load_metrics():
    path = os.path.join('results', 'model', 'metrics.json')
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except Exception:
        return None

st.title("Cardiac Variant Predictor (Demo)")
metrics = load_metrics()
if metrics:
    mcol1, mcol2, mcol3 = st.columns(3)
    mcol1.metric("ROC AUC", f"{metrics.get('ROC_AUC', 0):.3f}")
    mcol2.metric("PR AUC", f"{metrics.get('PR_AUC', 0):.3f}")
    mcol3.metric("Model", metrics.get('model_type', 'unknown'))

@st.cache_resource
def load_model():
    path = os.path.join('results', 'model', 'model_xgb.pkl')
    if not os.path.exists(path):
        st.warning("Model not found; please train first.")
        return None
    obj = joblib.load(path)
    return obj['model'], obj['features']

@st.cache_resource
def load_features():
    path = os.path.join('results', 'variant_features.csv')
    if not os.path.exists(path):
        st.warning("Features file not found; please run feature build.")
        return None
    return pd.read_csv(path)

model_features = load_model()
features_df = load_features()


def _rerun():
    """Streamlit rerun compatible with newer/older versions."""
    try:
        # Newer versions
        getattr(st, 'rerun')()
    except Exception:
        # Older versions
        if hasattr(st, 'experimental_rerun'):
            st.experimental_rerun()

col1, col2 = st.columns([2,1])
with col1:
    st.subheader("Lookup or Upload")
    choice = st.radio("Input method", ["Lookup by chrom:pos:ref:alt", "Upload small VCF (not implemented)"])

    chrom = st.text_input("Chromosome", value=str(st.session_state.get('chrom', '11')))
    pos = st.number_input("Position", value=int(st.session_state.get('pos', 2608632)), step=1)
    ref = st.text_input("Ref", value=str(st.session_state.get('ref', 'C')))
    alt = st.text_input("Alt", value=str(st.session_state.get('alt', 'CA')))

    if st.button("Predict"):
        if features_df is None or model_features is None:
            st.stop()
        model, feats = model_features
        q = features_df[(features_df['chrom'].astype(str)==str(chrom)) & (features_df['pos']==int(pos)) & (features_df['ref']==ref) & (features_df['alt']==alt)]
        if q.empty:
            st.error("Variant not found in precomputed features.")
        else:
            # build input row
            df = q.copy()
            mapping = {'D':1, 'P':1, 'A':1, 'T':0, 'B':0, 'N':0}
            for c in ['SIFT_pred','Polyphen2_HVAR_pred','Polyphen2_HDIV_pred']:
                if c in df.columns:
                    df[c] = df[c].map(mapping).fillna(0)
            cols = [c for c in ['gnomADg_AF','heart_lv_tpm','heart_aa_tpm','GERP++_RS','phyloP100way_vertebrate','SIFT_pred','Polyphen2_HVAR_pred','Polyphen2_HDIV_pred'] if c in df.columns]
            X = df[cols].fillna(0.0).values
            proba = float(model.predict_proba(X)[:,1][0])
            st.success(f"Predicted pathogenicity probability: {proba*100:.1f}%")

            # simple explanation using feature contributions from tree SHAP if available
            try:
                import shap
                explainer = shap.TreeExplainer(model)
                shap_vals = explainer.shap_values(X)

                # Select row and class robustly
                import numpy as np
                def _first_row(vals):
                    arr = vals
                    if isinstance(arr, list) and len(arr) > 0:
                        arr = arr[-1]  # assume last class is positive
                    if hasattr(arr, 'toarray'):
                        arr = arr.toarray()
                    arr = np.asarray(arr)
                    if arr.ndim == 2:
                        return arr[0]
                    return arr

                row_vals = _first_row(shap_vals)
                s = pd.Series(row_vals, index=cols).sort_values(key=abs, ascending=False)
                st.write("Top contributing features:")
                st.bar_chart(s.head(10))
                # brief evidence text
                evidence = []
                if 'gnomADg_AF' in df.columns:
                    af = float(df['gnomADg_AF'].iloc[0] or 0.0)
                    if af < 0.001:
                        evidence.append("Very low population AF")
                    elif af < 0.01:
                        evidence.append("Low population AF")
                if 'phyloP100way_vertebrate' in df.columns:
                    p = float(df['phyloP100way_vertebrate'].iloc[0] or 0.0)
                    if p > 2:
                        evidence.append("High phyloP conservation")
                if 'GERP++_RS' in df.columns:
                    g = float(df['GERP++_RS'].iloc[0] or 0.0)
                    if g > 2:
                        evidence.append("High GERP++ constraint")
                if 'SIFT_pred' in df.columns and df['SIFT_pred'].iloc[0] in ['D']:
                    evidence.append("SIFT damaging")
                if 'Polyphen2_HVAR_pred' in df.columns and df['Polyphen2_HVAR_pred'].iloc[0] in ['D','P']:
                    evidence.append("PolyPhen (HVAR) possibly/probably damaging")
                if evidence:
                    st.caption("Evidence: " + "; ".join(evidence))

                # SHAP CSV download (single variant)
                try:
                    import io
                    shap_df = pd.DataFrame({'feature': cols, 'shap_value': row_vals})
                    buf = io.StringIO()
                    shap_df.to_csv(buf, index=False)
                    st.download_button(
                        "Download SHAP CSV",
                        data=buf.getvalue(),
                        file_name=f"shap_{chrom}_{pos}_{ref}_{alt}.csv",
                        mime="text/csv",
                    )
                except Exception:
                    pass
            except Exception as e:
                st.info("SHAP explanation not available: "+str(e))

            # JSON download
            payload = {
                'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt,
                'prob_pathogenic': proba,
                'features_used': cols,
            }
            st.download_button(
                "Download prediction JSON",
                data=json.dumps(payload, indent=2),
                file_name=f"prediction_{chrom}_{pos}_{ref}_{alt}.json",
                mime="application/json",
            )

with col2:
    st.subheader("Demos")
    st.caption("Try examples from the precomputed feature table")
    demos = []
    try:
        # pick first pathogenic, first benign, and a random row
        if features_df is not None:
            dset = pd.read_csv(os.path.join('results','variant_dataset.csv'))
            pth = dset[dset['label']==1].head(1)
            ben = dset[dset['label']==0].head(1)
            vus = dset[dset['label'].isna()].head(1)
            for name, df0 in [("Pathogenic-like", pth), ("Benign-like", ben), ("VUS-like", vus)]:
                if not df0.empty:
                    row = df0.iloc[0]
                    demos.append((name, str(row['chrom']), int(row['pos']), str(row['ref']), str(row['alt'])))
    except Exception:
        pass
    for i, (label, c, p, r, a) in enumerate(demos):
        if st.button(label, key=f"demo_{i}"):
            st.session_state['chrom'] = c
            st.session_state['pos'] = p
            st.session_state['ref'] = r
            st.session_state['alt'] = a
            _rerun()
