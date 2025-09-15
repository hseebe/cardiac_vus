SHELL := /bin/zsh
.ONESHELL:

# Directories
ROOT := $(shell pwd)
DATA := $(ROOT)/data
SCRIPTS := $(ROOT)/scripts
RESULTS := $(ROOT)/results
BED := $(DATA)/bed
CLINVAR := $(ROOT)
DBNSFP := $(DATA)/dbnsfp
GTEX := $(DATA)/gtex
CONS := $(DATA)/conservation
VENV := $(ROOT)/.venv
PY := $(VENV)/bin/python
PIP := $(VENV)/bin/pip

# Files
CLINVAR_VCF := $(CLINVAR)/clinvar.vcf.gz
CLINVAR_TBI := $(CLINVAR)/clinvar.vcf.gz.tbi
CARDIAC_BED := $(BED)/cardiac_genes.bed
CARDIAC_VCF := $(ROOT)/cardiac_genes.vcf.gz
CARDIAC_CODING_VCF := $(ROOT)/cardiac_coding_variants.vcf.gz
ANNOTATED_VCF := $(RESULTS)/annotated_variants.vep.vcf
GTEX_GCT := $(GTEX)/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
GTEX_ATTR := $(GTEX)/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
GTEX_SUBSET := $(GTEX)/gtex_heart_subset.csv
PHYLOP := $(CONS)/hg38.phyloP100way.bw
PHASTCONS := $(CONS)/hg38.phastCons100way.bw
DBNSFP_ZIP := $(DBNSFP)/dbNSFP4.4a.zip
DBNSFP_MERGED := $(DBNSFP)/dbNSFP4.4a.txt.gz
DBNSFP_TBI := $(DBNSFP)/dbNSFP4.4a.txt.gz.tbi

# Default
.PHONY: all
all: setup download bed filter annotate gtex verify


.PHONY: setup
setup: $(VENV) deps
	@echo "Setup complete"

$(VENV):
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

.PHONY: deps
deps:
	@chmod +x scripts/*.sh || true
	@./scripts/check_deps.sh

.PHONY: vep
vep:
	@./scripts/ensure_vep.sh

.PHONY: download
download: download-clinvar download-dbnsfp download-gtex download-conservation

.PHONY: download-clinvar
download-clinvar: | $(DATA)
	@./scripts/download_clinvar.sh "$(CLINVAR_VCF)"

.PHONY: download-dbnsfp
download-dbnsfp: | $(DBNSFP)
	@./scripts/download_dbnsfp.sh "$(DBNSFP)"

.PHONY: download-gtex
download-gtex: | $(GTEX)
	@./scripts/download_gtex.sh "$(GTEX_GCT)"

.PHONY: download-conservation
download-conservation: | $(CONS)
	@./scripts/download_conservation.sh "$(CONS)"

.PHONY: bed
bed: | $(BED)
	$(PY) scripts/build_bed.py --out $(CARDIAC_BED)

.PHONY: filter
filter: $(CARDIAC_VCF) $(CARDIAC_CODING_VCF)

$(CARDIAC_VCF): $(CLINVAR_VCF) $(CARDIAC_BED)
	bcftools view -R $(CARDIAC_BED) $(CLINVAR_VCF) -Oz -o $@
	bcftools index -f $@

$(CARDIAC_CODING_VCF): $(CARDIAC_VCF)
	# Filter to missense/nonsense/frameshift if ANN or VEP fields exist; otherwise pass-through
	if bcftools view -h $(CARDIAC_VCF) | grep -E "ANN=|CSQ=" >/dev/null; then \
	  bcftools view -i 'INFO/ANN~"missense_variant|stop_gained|frameshift_variant|stop_lost" || INFO/CSQ~"missense_variant|stop_gained|frameshift_variant|stop_lost"' $(CARDIAC_VCF) -Oz -o $@; \
	  bcftools index -f $@; \
	else \
	  cp $(CARDIAC_VCF) $@ && bcftools index -f $@; \
	fi

.PHONY: annotate
annotate: | $(RESULTS)
	@./scripts/run_vep.sh "$(CARDIAC_CODING_VCF)" "$(ANNOTATED_VCF)" "$(DBNSFP_MERGED)"

.PHONY: gtex
gtex: $(GTEX_SUBSET)

$(GTEX_SUBSET): $(GTEX_GCT) $(VENV)
	ATTR_OPT=
	if [ -f "$(GTEX_ATTR)" ]; then ATTR_OPT="--attrs $(GTEX_ATTR)"; fi
	$(PY) scripts/extract_gtex_heart.py --gct $(GTEX_GCT) $$ATTR_OPT --out $(GTEX_SUBSET)

.PHONY: verify
verify:
	$(PY) scripts/verify_outputs.py \
	  --clinvar $(CLINVAR_VCF) \
	  --cardiac $(CARDIAC_VCF) \
	  --vep $(ANNOTATED_VCF) \
	  --gtex $(GTEX_SUBSET) \
	  --cons-phyloP $(PHYLOP) \
	  --cons-phastCons $(PHASTCONS) \
	  --dbnsfp-dir $(DBNSFP)

# Helpers
$(DATA):
	mkdir -p $(DATA)
$(RESULTS):
	mkdir -p $(RESULTS)
$(BED):
	mkdir -p $(BED)
$(DBNSFP):
	mkdir -p $(DBNSFP)
$(GTEX):
	mkdir -p $(GTEX)
$(CONS):
	mkdir -p $(CONS)

.PHONY: vep-cache
vep-cache:
	@./scripts/install_vep_cache.sh

# Clean up heavy build inputs but keep minimal demo runtime artifacts
.PHONY: clean-inference
clean-inference:
	@echo "Keeping only demo runtime artifacts..."
	@mkdir -p archive
	@if [ -d ".vep_cache" ]; then mv .vep_cache archive/; fi
	@if [ -d "data/conservation" ]; then mv data/conservation archive/; fi
	@if [ -f "data/dbnsfp/dbNSFP4.4a.zip" ]; then mv data/dbnsfp/dbNSFP4.4a.zip archive/; fi
	@find data/dbnsfp -type f -name 'dbNSFP4.4a_variant.chr*' -maxdepth 1 -exec mv {} archive/ \; || true
	@if [ -f "clinvar.vcf.gz" ]; then mv clinvar.vcf.gz* archive/; fi
	@if [ -f "cardiac_genes.vcf.gz" ]; then mv cardiac_genes.vcf.gz* archive/; fi
	@if [ -f "cardiac_coding_variants.vcf.gz" ]; then mv cardiac_coding_variants.vcf.gz* archive/; fi
	@if [ -f "results/annotated_variants.vep.vcf" ]; then rm -f results/annotated_variants.vep.vcf; fi
	@if [ -f "results/annotated_variants.vep.vcf_warnings.txt" ]; then mv results/annotated_variants.vep.vcf_warnings.txt archive/; fi
	@if [ -f "results/annotated_variants.vep.vcf_summary.html" ]; then mv results/annotated_variants.vep.vcf_summary.html archive/; fi
	@if [ -f "results/model/shap_values.csv" ]; then mv results/model/shap_values.csv archive/; fi
	@echo "Done. Minimal demo artifacts are preserved."

.PHONY: bundle
bundle:
	@echo "Creating slim deliverable tarball..."
	@bgzip -f results/annotated_variants.vep.vcf || true
	@tabix -p vcf results/annotated_variants.vep.vcf.gz || true
	tar -czf cardiac_vus_demo.tgz \
	  app.py README.md requirements.txt \
	  results/model \
	  results/variant_features.csv \
	  results/variant_dataset.csv \
	  results/annotated_variants.vep.vcf.gz \
	  results/annotated_variants.vep.vcf.gz.tbi \
	  data/gtex/gtex_heart_subset.csv
	@echo "Bundle: cardiac_vus_demo.tgz"
