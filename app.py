import streamlit as st
import pandas as pd
import joblib
import json
import os
import io
import time
import uuid
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Optional
import numpy as np

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


# ------------------------------
# Utilities for VCF upload flow
# ------------------------------

def _write_uploaded_file(uploaded_file, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    name = uploaded_file.name
    # Normalize extensions for safety
    if name.endswith('.vcf'):
        out_path = out_dir / 'input.vcf'
        with open(out_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        return out_path
    elif name.endswith('.vcf.gz'):
        gz_path = out_dir / 'input.vcf.gz'
        with open(gz_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        # Decompress to .vcf (VEP can read gz, but we normalize here for consistency)
        out_path = out_dir / 'input.vcf'
        with gzip.open(gz_path, 'rb') as fin, open(out_path, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        return out_path
    else:
        # Attempt best-effort treat as plain text VCF
        out_path = out_dir / 'input.vcf'
        with open(out_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())
        return out_path


def _run_vep(in_vcf: Path, out_vcf: Path) -> tuple[bool, str]:
    """Run VEP via provided helper script with Docker fallback. Returns (ok, log)."""
    script = Path('scripts/run_vep.sh')
    dbnsfp = Path('data/dbnsfp/dbNSFP4.4a.txt.gz')
    cmd = ['zsh', str(script), str(in_vcf), str(out_vcf)]
    if dbnsfp.exists():
        cmd.append(str(dbnsfp))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        ok = proc.returncode == 0
        log = (proc.stdout or '') + '\n' + (proc.stderr or '')
        return ok, log
    except Exception as e:
        return False, f"Failed to invoke VEP: {e}"


def _parse_csq_header(vcf_path: Path):
    import re
    fields = []
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    with opener(vcf_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                m = re.search(r'Format: (.*)\">', line)
                if m:
                    fields = m.group(1).split('|')
            if line.startswith('#CHROM'):
                break
    return fields


def _build_features_from_vep(vcf_path: Path, gtex_csv: Path, clinvar_vcf: Optional[Path] = None) -> pd.DataFrame:
    """Lightweight feature builder compatible with scripts/build_features.py but clinvar is optional."""
    import re

    csq_fields = _parse_csq_header(vcf_path)
    want = [
        'Consequence', 'IMPACT', 'SYMBOL', 'gnomADg_AF',
        'SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred',
        'GERP++_RS', 'phyloP100way_vertebrate',
    ]

    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    records = []
    with opener(vcf_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt = parts[:5]
            info = parts[7]
            m = re.search(r'CSQ=([^;]+)', info)
            if not m:
                continue
            best = None
            for csq in m.group(1).split(','):
                vals = csq.split('|')
                if len(vals) != len(csq_fields):
                    continue
                row = {k: v for k, v in zip(csq_fields, vals)}
                if row.get('CANONICAL') == 'YES':
                    best = row
                    break
                if best is None:
                    best = row
            if best is None:
                continue
            rec = {
                'chrom': chrom.lstrip('chr'),
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'id': vid,
            }
            for k in want:
                rec[k] = best.get(k)
            records.append(rec)

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return df
    for col in ['gnomADg_AF', 'GERP++_RS', 'phyloP100way_vertebrate']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Join GTEx by gene symbol
    try:
        gtex = pd.read_csv(gtex_csv)
        rename_map = {}
        for c in gtex.columns:
            lc = c.strip().lower()
            if lc.startswith('heart - left ventricle'):
                rename_map[c] = 'heart_lv_tpm'
            if lc.startswith('heart - atrial appendage'):
                rename_map[c] = 'heart_aa_tpm'
        if 'gene' not in gtex.columns:
            if 'Description' in gtex.columns:
                gtex = gtex.rename(columns={'Description': 'gene'})
            elif 'Name' in gtex.columns:
                gtex = gtex.rename(columns={'Name': 'gene'})
        gtex = gtex.rename(columns=rename_map)
        cols = ['gene'] + [c for c in ['heart_lv_tpm', 'heart_aa_tpm'] if c in gtex.columns]
        gtex = gtex[cols].drop_duplicates('gene')
        df = df.merge(gtex, left_on='SYMBOL', right_on='gene', how='left').drop(columns=['gene'], errors='ignore')
    except Exception:
        pass

    # Optional ClinVar labels if available
    if clinvar_vcf is not None and Path(clinvar_vcf).exists():
        labels = {}
        opener2 = gzip.open if str(clinvar_vcf).endswith('.gz') else open
        with opener2(clinvar_vcf, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                chrom, pos, _vid, ref, alt = parts[:5]
                info = parts[7]
                m = re.search(r'CLNSIG=([^;]+)', info)
                if not m:
                    continue
                sigs = {s.strip() for s in m.group(1).split(',')}
                label = None
                if 'Pathogenic' in sigs or 'Likely_pathogenic' in sigs:
                    label = 1
                elif 'Benign' in sigs or 'Likely_benign' in sigs:
                    label = 0
                if label is not None:
                    labels[(chrom.lstrip('chr'), int(pos), ref, alt)] = label
        df['label'] = [labels.get((str(ch), int(p), r, a)) if pd.notnull(p) else None for ch, p, r, a in zip(df['chrom'], df['pos'], df['ref'], df['alt'])]

    return df


def _predict_on_df(df: pd.DataFrame, model, feature_cols: list[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    mapping = {'D': 1, 'P': 1, 'A': 1, 'T': 0, 'B': 0, 'N': 0}
    for c in ['SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred']:
        if c in df.columns:
            df[c] = df[c].map(mapping).fillna(0)
    cols = [c for c in ['gnomADg_AF', 'heart_lv_tpm', 'heart_aa_tpm', 'GERP++_RS', 'phyloP100way_vertebrate', 'SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred'] if c in df.columns]
    if not cols:
        return df.assign(prob_pathogenic=np.nan), np.empty((0,)), []
    X = df[cols].fillna(0.0).values
    proba = model.predict_proba(X)[:, 1]
    df = df.copy()
    df['prob_pathogenic'] = proba
    return df, proba, cols


def _evidence_bullets(row: pd.Series) -> list[str]:
    ev = []
    try:
        af = float(row.get('gnomADg_AF') or 0.0)
        if af == 0 or af < 1e-5:
            ev.append('PM2_supporting: absent/ultra-rare in gnomAD')
        elif af < 0.001:
            ev.append('PM2_moderate: very low population AF')
        elif af < 0.01:
            ev.append('BS1_supporting: low AF')
    except Exception:
        pass
    try:
        p = float(row.get('phyloP100way_vertebrate') or 0.0)
        if p > 2:
            ev.append('Conserved site (phyloP high)')
    except Exception:
        pass
    try:
        g = float(row.get('GERP++_RS') or 0.0)
        if g > 2:
            ev.append('Constraint (GERP++ high)')
    except Exception:
        pass
    if str(row.get('SIFT_pred', '')).upper() == 'D':
        ev.append('PP3: SIFT damaging')
    if str(row.get('Polyphen2_HVAR_pred', '')).upper() in {'D', 'P'}:
        ev.append('PP3: PolyPhen HVAR possibly/probably damaging')
    if row.get('IMPACT') in {'HIGH'}:
        ev.append('PVS1_suggestive: high-impact consequence')
    return ev

col1, col2 = st.columns([2, 1])
with col1:
    st.subheader("Lookup or Upload")
    choice = st.radio("Input method", ["Lookup by chrom:pos:ref:alt", "Upload small VCF"])

    chrom = st.text_input("Chromosome", value=str(st.session_state.get('chrom', '11')))
    pos = st.number_input("Position", value=int(st.session_state.get('pos', 2608632)), step=1)
    ref = st.text_input("Ref", value=str(st.session_state.get('ref', 'C')))
    alt = st.text_input("Alt", value=str(st.session_state.get('alt', 'CA')))

    if choice == "Lookup by chrom:pos:ref:alt" and st.button("Predict"):
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

    if choice == "Upload small VCF":
        st.markdown("Upload a small VCF (dozens–hundreds of variants). We'll annotate with VEP, build features, and predict.")
        up = st.file_uploader("VCF file", type=["vcf", "vcf.gz"], accept_multiple_files=False)
        run = st.button("Annotate and Predict")
        demo = st.button("Run example VCF (data/dbnsfp/try.vcf)")
        def _process_vcf(file_like_path: Path):
            if model_features is None:
                st.warning("Model not found; please train first.")
                st.stop()
            model, feats = model_features
            job_id = time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:6]
            workdir = Path('results') / 'uploads' / job_id
            with st.status("Processing VCF…", expanded=True) as status:
                st.write("Preparing input…")
                workdir.mkdir(parents=True, exist_ok=True)
                if isinstance(file_like_path, Path) and file_like_path.exists():
                    in_vcf = workdir / 'input.vcf'
                    if str(file_like_path).endswith('.gz'):
                        with gzip.open(file_like_path, 'rb') as fin, open(in_vcf, 'wb') as fout:
                            shutil.copyfileobj(fin, fout)
                    else:
                        shutil.copy(file_like_path, in_vcf)
                else:
                    st.error("Input VCF not found.")
                    status.update(label="Failed", state="error")
                    st.stop()
                st.write("Annotating with VEP…")
                out_vcf = workdir / 'annotated.vep.vcf'
                ok, log = _run_vep(in_vcf, out_vcf)
                if not ok:
                    st.error("VEP annotation failed. Install VEP locally (scripts/ensure_vep.sh) or ensure Docker is available. Logs below.")
                    st.code(log)
                    status.update(label="Failed", state="error")
                    st.stop()
                st.write("Building features…")
                gtex_csv = Path('data/gtex/gtex_heart_subset.csv')
                clinvar_vcf = Path('cardiac_genes.vcf.gz') if Path('cardiac_genes.vcf.gz').exists() else None
                df = _build_features_from_vep(out_vcf, gtex_csv, clinvar_vcf)
                if df.empty:
                    st.error("No variants with usable annotations were found.")
                    status.update(label="No usable variants", state="error")
                    st.stop()
                st.write(f"Feature rows: {len(df)}")
                pred_df, proba, used_cols = _predict_on_df(df, model, feats)
                if not used_cols:
                    st.warning("Annotated VCF lacks required features (dbNSFP/gnomAD). Try re-running with dbNSFP plugin available. See scripts/run_vep.sh.")
                pred_path = workdir / 'predictions.csv'
                pred_df.to_csv(pred_path, index=False)
                status.update(label="Done", state="complete")
            return workdir, out_vcf, pred_df, used_cols

        if run:
            if up is None:
                st.warning("Please choose a VCF file first.")
                st.stop()
            # Save and process uploaded file
            job_id = time.strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:6]
            workdir = Path('results') / 'uploads' / job_id
            in_vcf = _write_uploaded_file(up, workdir)
            workdir, out_vcf, pred_df, used_cols = _process_vcf(in_vcf)

            # Display results
            st.success(f"Annotated and scored {len(pred_df)} variants. Used features: {', '.join(used_cols) if used_cols else 'none'}")

            # Evidence bullets per variant
            pred_df['evidence'] = pred_df.apply(_evidence_bullets, axis=1)
            show_cols = ['chrom', 'pos', 'ref', 'alt', 'SYMBOL', 'Consequence', 'IMPACT', 'gnomADg_AF', 'heart_lv_tpm', 'heart_aa_tpm', 'GERP++_RS', 'phyloP100way_vertebrate', 'prob_pathogenic', 'evidence']
            show = [c for c in show_cols if c in pred_df.columns]
            st.dataframe(pred_df[show].sort_values('prob_pathogenic', ascending=False), use_container_width=True)

            # Downloads
            with open(out_vcf, 'r', encoding='utf-8', errors='ignore') as f:
                ann_text = f.read()
            st.download_button("Download annotated VCF", data=ann_text, file_name=f"annotated_{job_id}.vcf", mime="text/vcf")
            st.download_button("Download predictions CSV", data=pred_df.to_csv(index=False), file_name=f"predictions_{job_id}.csv", mime="text/csv")

            # Scientific wow: SHAP and interactive views
            try:
                import shap
                if used_cols:
                    X = pred_df[used_cols].fillna(0.0).values
                    explainer = shap.TreeExplainer(model)
                    vals = explainer.shap_values(X)
                    # For binary models, take positive class
                    if isinstance(vals, list) and len(vals) > 1:
                        vals = vals[-1]
                    shap_abs_mean = np.abs(np.asarray(vals)).mean(axis=0)
                    feat_imp = pd.Series(shap_abs_mean, index=used_cols).sort_values(ascending=False)
                    st.subheader("Global feature importance (SHAP)")
                    st.bar_chart(feat_imp.head(15))
                    # Per-variant top drivers (first 1-3 variants)
                    st.caption("Per-variant drivers (first up to 3 variants)")
                    n_show = min(3, len(pred_df))
                    for i in range(n_show):
                        row_vals = np.asarray(vals)[i]
                        s = pd.Series(row_vals, index=used_cols).sort_values(key=np.abs, ascending=False).head(8)
                        st.write(f"Variant {i+1}: {pred_df.loc[pred_df.index[i], 'chrom']}:{pred_df.loc[pred_df.index[i], 'pos']} {pred_df.loc[pred_df.index[i], 'ref']}>{pred_df.loc[pred_df.index[i], 'alt']}  — P(path)={pred_df.loc[pred_df.index[i], 'prob_pathogenic']:.2f}")
                        st.bar_chart(s)
            except Exception as e:
                st.info(f"SHAP summary unavailable: {e}")

            # Interactive scatter: probability vs conservation
            try:
                import altair as alt
                base = pred_df.copy()
                base['-log10_AF'] = -np.log10(base['gnomADg_AF'].replace(0, np.nan)) if 'gnomADg_AF' in base.columns else np.nan
                xfield = 'phyloP100way_vertebrate' if 'phyloP100way_vertebrate' in base.columns else ('GERP++_RS' if 'GERP++_RS' in base.columns else None)
                if xfield is not None:
                    chart = alt.Chart(base).mark_circle(opacity=0.8).encode(
                        x=alt.X(xfield, title=f"{xfield} (conservation)"),
                        y=alt.Y('prob_pathogenic', title='Predicted pathogenicity'),
                        color=alt.Color('Consequence:N', legend=None),
                        tooltip=['chrom', 'pos', 'ref', 'alt', 'SYMBOL', 'Consequence', 'IMPACT', 'prob_pathogenic']
                    ).interactive()
                    st.subheader("Probability vs conservation")
                    st.altair_chart(chart, use_container_width=True)
            except Exception:
                pass

        if demo:
            demo_path = Path('data/dbnsfp/try.vcf')
            if not demo_path.exists():
                st.error("Demo VCF not found at data/dbnsfp/try.vcf")
                st.stop()
            workdir, out_vcf, pred_df, used_cols = _process_vcf(demo_path)
            st.success(f"Annotated and scored {len(pred_df)} variants from demo VCF.")
            show_cols = ['chrom', 'pos', 'ref', 'alt', 'SYMBOL', 'Consequence', 'IMPACT', 'gnomADg_AF', 'heart_lv_tpm', 'heart_aa_tpm', 'GERP++_RS', 'phyloP100way_vertebrate', 'prob_pathogenic']
            show = [c for c in show_cols if c in pred_df.columns]
            st.dataframe(pred_df[show].sort_values('prob_pathogenic', ascending=False), use_container_width=True)

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
