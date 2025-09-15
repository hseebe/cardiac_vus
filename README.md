# Cardiac Variant Explainable AI – Data Prep Pipeline

This workspace sets up a reproducible pipeline to:
- Ensure system and Python dependencies
- Download required datasets (ClinVar, dbNSFP, GTEx, conservation bigWigs)
- Build a gene BED for a cardiac gene panel
- Filter ClinVar to the gene panel and coding-impact variants
- Annotate variants with Ensembl VEP (offline) and dbNSFP plugin
- Extract heart tissue expression from GTEx
- Verify deliverables

## Gene panel (GRCh38)
- MYH7, MYBPC3, TNNT2, TTN, SCN5A, KCNQ1

## Quick start

1) Setup Python env and check CLI tools
- Creates .venv and installs Python deps (pandas, pyBigWig)
- Ensures wget/curl, bcftools, tabix, unzip; installs via Homebrew if available

2) Download datasets
- ClinVar VCF + index
- dbNSFP 4.4a and prepare a single tabix-indexed file for VEP plugin
- GTEx median TPM per tissue GCT
- UCSC 100-way phyloP and phastCons bigWigs (hg38)

3) Build BED for cardiac gene panel

4) Filter & annotate
- Filter ClinVar by BED
- Annotate with VEP (offline cache) + dbNSFP plugin

5) Extract GTEx subset (Heart: Left Ventricle, Atrial Appendage)

6) Verify deliverables

## Make targets
- make setup – Python venv + check/install CLI deps
- make download – Download datasets
- make bed – Build `data/bed/cardiac_genes.bed`
- make filter – Filter ClinVar to BED → `cardiac_genes.vcf.gz`
- make annotate – VEP + dbNSFP → `results/annotated_variants.vep.vcf`
- make gtex – Extract heart TPM subset → `data/gtex/gtex_heart_subset.csv`
- make verify – Check all expected files
- make all – End-to-end

## Notes
- VEP cache: The pipeline uses VEP offline cache for GRCh38. If not present, run `make vep-cache`.
- dbNSFP plugin: The pipeline prepares a bgzipped+tabix’d merged dbNSFP 4.4a file for the VEP `dbNSFP` plugin.
- GTEx median: Although the instructions mention per-sample GCT, we use the published median TPM per tissue file for direct extraction.
 - dbNSFP access: If the FTP login fails (credentials change periodically), download `dbNSFP4.4a.zip` manually from the official source into `data/dbnsfp/` and rerun `make download-dbnsfp`.

## Outputs
- Root: `clinvar.vcf.gz`, `clinvar.vcf.gz.tbi`, `cardiac_genes.vcf.gz`, `cardiac_genes.vcf.gz.tbi`
- results/: `annotated_variants.vep.vcf`
- data/gtex/: `gtex_heart_subset.csv`
- data/conservation/: `hg38.phyloP100way.bw`, `hg38.phastCons100way.bw`
- data/dbnsfp/: Prepared dbNSFP file(s) and index

## Modeling and demo

After `make annotate` and `make gtex` have produced the feature inputs, run:

1) Build labeled dataset
- Creates `results/variant_dataset.csv`

2) Train classifier and SHAP
- Creates `results/model/model_xgb.pkl` (or a fallback model on macOS), `results/model/metrics.json`, and SHAP files.

3) Launch the Streamlit demo
- Opens an interactive page to query a variant present in the features and view a probability and simple explanation.

Notes:
- On macOS without libomp, XGBoost may not load. The training script will automatically fall back to scikit-learn's GradientBoostingClassifier and still write `model_xgb.pkl` for the app to consume. To enable native XGBoost, install OpenMP via Homebrew: `brew install libomp`.

## Setup and run

1) Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2) Build features and dataset (if not already present)

```bash
# Requires prior VEP annotation and GTEx subset
.venv/bin/python scripts/build_dataset.py \
	--features results/variant_features.csv \
	--clinvar clinvar.vcf.gz \
	--out results/variant_dataset.csv

.venv/bin/python scripts/train_xgb.py \
	--data results/variant_dataset.csv \
	--outdir results/model
```

3) Launch the app

```bash
source .venv/bin/activate
streamlit run app.py
```

## Data sources and licenses

- ClinVar: NCBI ClinVar VCF (see NCBI terms of use)
- dbNSFP 4.4a: Database of Non-synonymous Functional Predictions (documentation and license in `data/dbnsfp/`)
- GTEx v8: The GTEx Consortium (see GTEx Portal; use per their terms)
- Conservation bigWigs: UCSC phyloP100way and phastCons100way (see UCSC terms)

## Known caveats

- dbNSFP subset: The dbNSFP file is a subset matching variants in our VCF to keep the plugin fast/light.
- macOS XGBoost: If `libomp` is not installed, training falls back to scikit-learn. Install with `brew install libomp` to enable native XGBoost.