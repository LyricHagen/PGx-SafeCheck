from typing import Dict, List, Any

# -----------------
# Parsing helpers
# -----------------

def parse_23andme_raw(text: str) -> Dict[str, str]:
    data = {}
    for line in text.splitlines():
        if not line or line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 4:
            parts = line.strip().split(',')
        if len(parts) < 4:
            continue
        rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
        data[rsid] = genotype.strip()
    return data

def parse_vcf(text: str) -> Dict[str, str]:
    rs_to_gt = {}
    sample_index = None
    for line in text.splitlines():
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_index = 9 if len(header) >= 10 else None
            continue
        if not line or line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        if rsid == '.' or not rsid.startswith('rs'):
            continue
        alleles = [ref] + alt.split(',')
        gt = None
        if sample_index is not None and len(parts) > sample_index:
            fmt = parts[8].split(':')
            sample = parts[sample_index].split(':')
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
                gt_val = sample[gt_idx] if gt_idx < len(sample) else None
                if gt_val:
                    sep = '/' if '/' in gt_val else '|' if '|' in gt_val else None
                    if sep:
                        a, b = gt_val.split(sep)
                        try:
                            gt = alleles[int(a)] + alleles[int(b)]
                        except Exception:
                            gt = None
        if gt is None:
            gt = alleles[0] + alleles[0]
        rs_to_gt[rsid] = gt
    return rs_to_gt

# -----------------
# PGx interpretation
# -----------------

def count_alt(geno: str, alt_allele: str) -> int:
    if not geno or geno == '--' or len(geno) < 2:
        return 0
    return int(geno[0] == alt_allele) + int(geno[1] == alt_allele)

def call_cyp2c19(rs: Dict[str, str]) -> Dict[str, Any]:
    g2 = rs.get('rs4244285')
    g3 = rs.get('rs4986893')
    g17 = rs.get('rs12248560')
    counts = {
        '*2': count_alt(g2 or '', 'A'),
        '*3': count_alt(g3 or '', 'A'),
        '*17': count_alt(g17 or '', 'T'),
    }
    alleles = []
    for star in ['*2','*3','*17']:
        alleles += [star] * min(counts[star], 2 - len(alleles))
        if len(alleles) >= 2:
            break
    while len(alleles) < 2:
        alleles.append('*1')

    diplotype = '/'.join(alleles)
    phenotype = 'Normal metabolizer'
    if alleles.count('*2') + alleles.count('*3') == 2:
        phenotype = 'Poor metabolizer'
    elif (('*2' in alleles) or ('*3' in alleles)) and ('*17' in alleles or '*1' in alleles):
        phenotype = 'Intermediate metabolizer'
    elif alleles == ['*1','*1']:
        phenotype = 'Normal metabolizer'
    elif alleles.count('*17') == 1 and '*1' in alleles:
        phenotype = 'Rapid metabolizer'
    elif alleles.count('*17') == 2:
        phenotype = 'Ultrarapid metabolizer'

    caveats = [
        "Approximate call using only *2, *3, *17. Other alleles and phasing not assessed.",
        "Phenotype mapping simplified for demo; consult CPIC for authoritative guidance.",
    ]
    associations = [
        {
            "area": "Antiplatelet therapy",
            "drug_examples": ["clopidogrel"],
            "summary": (
                "CYP2C19 {} may have altered activation of clopidogrel, potentially changing antiplatelet effect."
                .format(phenotype.lower())
            ),
            "strength": "well-established",
            "notes": [
                "Discuss antiplatelet strategy with a clinician; do not change therapy based on this demo.",
            ],
        }
    ]
    return {
        "gene": "CYP2C19",
        "diplotype": diplotype,
        "phenotype": f"{phenotype} (approximate)",
        "caveats": caveats,
        "associations": associations,
    }

def call_slco1b1(rs: Dict[str, str]) -> Dict[str, Any]:
    g = rs.get('rs4149056', 'TT')
    c_count = count_alt(g, 'C')
    category = "Typical function"
    if c_count == 1:
        category = "Decreased function (heterozygous)"
    elif c_count == 2:
        category = "Decreased function (homozygous)"
    caveats = [
        "Single-variant approximation; full haplotypes and other variants not assessed.",
        "Risk depends on statin type/dose and patient-specific factors.",
    ]
    associations = [
        {
            "area": "Statins (e.g., simvastatin)",
            "drug_examples": ["simvastatin"],
            "summary": "The SLCO1B1 521C variant is associated with higher statin plasma levels and increased myopathy risk.",
            "strength": "well-established",
            "notes": [
                "Discuss statin choice/dose with a clinician; do not start/stop medication based on this demo."
            ],
        }
    ]
    return {
        "gene": "SLCO1B1",
        "genotype": g,
        "functional_category": category + " (approximate)",
        "caveats": caveats,
        "associations": associations,
    }

def call_hlab_5701_proxy(rs: Dict[str, str]) -> Dict[str, Any]:
    g = rs.get('rs2395029', '--')
    g_count = count_alt(g, 'G')
    tag_status = "Proxy not present"
    if g_count == 1:
        tag_status = "Possible proxy (heterozygous)"
    elif g_count == 2:
        tag_status = "Possible proxy (homozygous)"
    caveats = [
        "rs2395029 is an imperfect proxy for HLA-B*57:01; confirmation testing is required for clinical use.",
        "Tag performance varies by ancestry; do not make treatment decisions based on this demo.",
    ]
    associations = [
        {
            "area": "Abacavir hypersensitivity",
            "drug_examples": ["abacavir"],
            "summary": "Presence of proxy allele indicates potential abacavir hypersensitivity risk.",
            "strength": "well-established (for the allele); proxy tag is imperfect",
            "notes": ["Discuss with a clinician and consider confirmatory HLA-B*57:01 testing if relevant."],
        }
    ]
    return {
        "gene": "HLA-B*57:01 (proxy via rs2395029)",
        "proxy_status": tag_status,
        "genotype": g,
        "caveats": caveats,
        "associations": associations,
    }

# -----------------
# Orchestrator
# -----------------

KNOWN_RSIDS = {
    'rs4244285',  # CYP2C19 *2
    'rs4986893',  # CYP2C19 *3
    'rs12248560', # CYP2C19 *17
    'rs4149056',  # SLCO1B1
    'rs2395029',  # HLA-B*57:01 proxy
}

def extract_known_rsids(rs_map: Dict[str, str]) -> Dict[str, str]:
    return {k: v for k, v in rs_map.items() if k in KNOWN_RSIDS}

def analyze(rs_map: Dict[str, str]) -> Dict[str, Any]:
    subset = extract_known_rsids(rs_map)
    results: List[Dict[str, Any]] = []
    results.append(call_cyp2c19(subset))
    results.append(call_slco1b1(subset))
    results.append(call_hlab_5701_proxy(subset))

    summary_points = [
        "This report is informational only and not a substitute for professional medical advice.",
        "Associations are approximate; confirm with validated clinical testing and consult guidelines.",
        "Do not change medications based on this report.",
    ]
    return {
        "summary": summary_points,
        "results": results,
        "observed_rsids": subset,
        "missing_rsids": sorted(list(KNOWN_RSIDS - set(subset.keys()))),
        "version": "demo-0.1.0",
    }
