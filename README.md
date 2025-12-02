# PGx SafeCheck

Minimal **pharmacogenomics (PGx)** app: upload genetic data, get **informative, non-prescriptive** gene–drug associations.  
Not a medical device; strictly **educational/decision-support**.

**⚠️ Safety**
- **Not medical advice** — do not change medications based on this app.
- Findings are **limited**; CNVs, rare variants, phased haplotypes, and ancestry-specific effects are not fully handled.
- Focused on **CYP2C19, SLCO1B1, HLA-B*57:01 (proxy)**. CYP2D6 is omitted (complex).

## Quickstart

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
uvicorn backend.app:app --reload --port 8000
open http://localhost:8000

# Features

- Upload 23andMe raw text or simple VCF with rsIDs.

**Parses:**
- **CYP2C19:** rs4244285 (*2), rs4986893 (*3), rs12248560 (*17)  
- **SLCO1B1:** rs4149056  
- **HLA-B*57:01 (proxy):** rs2395029  

**Returns:**
- Approx diplotype/phenotype  
- Plain-language interpretation  
- Alerts for potential altered drug response  
- Caveats about limitations  

**Example:**  
Upload 23andMe file → gets:  
- **CYP2C19:** *1/*2, Intermediate metabolizer → reduced clopidogrel activation  
- **SLCO1B1:** TC, decreased function → higher statin myopathy risk  
- **HLA-B*57:01:** proxy absent → limited info  

# Limitations

- No prescribing or dosing guidance  
- Not guaranteed accurate for all populations  
- Informational only  

# Sources

- CPIC  
- PharmGKB  
- FDA PGx Table & labeling (DailyMed/OpenFDA)  

# Sample files

- `sample_data/23andme_example.txt`  
- `sample_data/simple.vcf`  

# Extending

- Add genes in `backend/pgx_logic.py`; keep outputs non-prescriptive.  
- Clinical use requires validation, QMS, risk management (ISO 14971), oversight, licensed knowledge base.
