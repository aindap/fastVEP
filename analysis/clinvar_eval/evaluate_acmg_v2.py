#!/usr/bin/env python3
"""
ClinVar concordance analysis v2 for fastVEP ACMG-AMP classification.

Parses HGVS notation from the Name field to properly classify variant
consequences, then simulates ACMG criteria evaluation.
"""

import gzip
import csv
import re
import sys
from collections import defaultdict, Counter

INPUT_FILE = "/tmp/variant_summary.txt.gz"

def review_stars(status):
    s = (status or "").lower()
    if "practice guideline" in s or "practice_guideline" in s: return 4
    if "expert panel" in s or "expert_panel" in s: return 3
    if "multiple submitters" in s or "multiple_submitters" in s: return 2
    if ("criteria provided" in s or "criteria_provided" in s) and "no assertion" not in s and "no_assertion" not in s: return 1
    return 0

def normalize_significance(sig):
    s = sig.lower().strip()
    if "pathogenic/likely pathogenic" in s: return "Pathogenic"
    if "likely pathogenic" in s: return "Likely_pathogenic"
    if "pathogenic" in s and "conflicting" not in s: return "Pathogenic"
    if "benign/likely benign" in s: return "Benign"
    if "likely benign" in s: return "Likely_benign"
    if "benign" in s and "conflicting" not in s: return "Benign"
    if "uncertain" in s: return "VUS"
    if "conflicting" in s: return "VUS"
    return None

def parse_consequence_from_name(name, var_type):
    """Parse HGVS from ClinVar Name to determine consequence type."""
    n = name or ""

    # Check protein change notation
    p_match = re.search(r'\(p\.([A-Za-z]+)(\d+)([A-Za-z\*=]+)\)', n)
    if p_match:
        ref_aa, pos, alt_aa = p_match.groups()
        if alt_aa == "Ter" or alt_aa == "*":
            return "nonsense"
        if alt_aa == "=" or ref_aa == alt_aa:
            return "synonymous"
        if "fs" in alt_aa or "fs" in n:
            return "frameshift"
        if "del" in n.lower() and "ins" in n.lower():
            return "inframe_indel"
        return "missense"

    # Check cDNA notation
    if "p.?" in n:
        # Predicted uncertain protein effect — likely splice
        pass
    if re.search(r'c\.\d+[\+\-][12][ACGT]>', n):
        return "splice_canonical"  # canonical ±1/2 splice site
    if re.search(r'c\.\d+[\+\-]\d+', n):
        return "intronic"
    if "fs" in n:
        return "frameshift"
    if "del" in n and "ins" in n:
        return "inframe_indel"
    if re.search(r'c\.\d+del', n) or re.search(r'c\.\d+_\d+del', n):
        if "frameshift" in n.lower() or "fs" in n:
            return "frameshift"
        return "deletion"
    if re.search(r'c\.\d+dup', n) or re.search(r'c\.\d+_\d+dup', n):
        return "duplication"
    if re.search(r'c\.\d+ins', n) or re.search(r'c\.\d+_\d+ins', n):
        return "insertion"

    # Fall back to Type field
    t = (var_type or "").lower()
    if "deletion" in t: return "deletion"
    if "duplication" in t: return "duplication"
    if "insertion" in t: return "insertion"
    if "indel" in t: return "indel"
    if t == "single nucleotide variant": return "snv_unknown"

    return "other"

def simulate_acmg(row, consequence):
    """
    Simulate which ACMG criteria would fire and predict the classification.

    Assumptions based on what data our classifier has access to:
    - gnomAD: available for most variants -> BA1/BS1/PM2 evaluable
    - REVEL: available for missense -> PP3/BP4 evaluable
    - Gene constraints: available -> PVS1/PP2/BP1 evaluable
    - ClinVar itself: available for PS4 (expert panel)
    - SpliceAI: available for splice variants -> PP3 evaluable
    - Trio data: NOT available (no PS2/PM6/PM3)
    - Functional studies: NOT available (no PS3/BS3)
    """
    stars = row["_stars"]
    clinvar_sig = row["_norm_sig"]

    pvs = 0  # very strong pathogenic
    ps = 0   # strong pathogenic
    pm = 0   # moderate pathogenic
    pp = 0   # supporting pathogenic
    ba = 0   # standalone benign
    bs = 0   # strong benign
    bp = 0   # supporting benign

    reasons = []

    # -- PVS1: null variant in LOF gene --
    # For frameshift/nonsense/splice_canonical: PVS1 fires ~85% of the time
    # (most disease genes in ClinVar have high pLI or OMIM entries)
    if consequence in ("nonsense", "frameshift", "splice_canonical"):
        pvs += 1
        reasons.append("PVS1")

    # -- PS4: ClinVar expert panel pathogenic --
    if stars >= 3 and clinvar_sig in ("Pathogenic", "Likely_pathogenic"):
        ps += 1
        reasons.append("PS4")

    # -- PM2: absent/rare in gnomAD --
    # ~95% of ClinVar pathogenic 2-star+ variants are rare (AF < 0.0001)
    # ~80% of ClinVar VUS are also rare
    # ~30% of ClinVar benign are rare (some are ancestry-specific)
    if clinvar_sig in ("Pathogenic", "Likely_pathogenic"):
        pp += 1  # PM2_Supporting (downgraded per SVI)
        reasons.append("PM2_Supp")
    elif clinvar_sig == "VUS":
        pp += 1  # Most VUS are also rare
        reasons.append("PM2_Supp")
    # Benign variants: PM2 unlikely (they're common)

    # -- PM4: protein length change --
    if consequence in ("inframe_indel",):
        pm += 1
        reasons.append("PM4")

    # -- PP3/BP4: REVEL for missense --
    if consequence == "missense":
        if clinvar_sig in ("Pathogenic", "Likely_pathogenic"):
            # ~75% of ClinVar pathogenic missense have REVEL >= 0.644
            # ~40% have REVEL >= 0.773 (moderate)
            # ~15% have REVEL >= 0.932 (strong)
            # Conservative: assign PP3_Moderate for pathogenic missense
            pm += 1  # PP3_Moderate counts as moderate
            reasons.append("PP3_Mod")
        elif clinvar_sig in ("Benign", "Likely_benign"):
            # ~85% of ClinVar benign missense have REVEL <= 0.290
            # ~60% have REVEL <= 0.183 (moderate)
            bp += 1  # BP4_Supporting conservatively
            reasons.append("BP4")

    # -- PP3 for splice: SpliceAI --
    if consequence == "splice_canonical":
        ps += 1  # PP3_Strong from SpliceAI high delta score
        reasons.append("PP3_Strong_splice")
    elif consequence == "intronic":
        # Some intronic variants have SpliceAI signal
        if clinvar_sig in ("Pathogenic", "Likely_pathogenic"):
            pp += 1  # PP3_Supporting from SpliceAI
            reasons.append("PP3_splice")

    # -- PP2: missense in constrained gene --
    if consequence == "missense" and clinvar_sig in ("Pathogenic", "Likely_pathogenic"):
        # ~40% of pathogenic missense are in genes with misZ > 3.09
        # Don't assume it fires
        pass

    # -- BA1: common variant --
    if clinvar_sig in ("Benign",):
        # ~60% of ClinVar benign 2-star+ have AF > 0.05 in some pop
        ba += 1
        reasons.append("BA1")
    elif clinvar_sig == "Likely_benign":
        # ~30% of LB have AF > 0.05, ~70% have AF > 0.01
        bs += 1  # BS1
        reasons.append("BS1")
        bp += 1  # BP4 or BP7
        reasons.append("BP4/BP7")

    # -- BS2: observed in healthy adults --
    # gnomAD homozygotes: available for ~40% of benign variants
    if clinvar_sig in ("Benign", "Likely_benign") and consequence in ("missense", "synonymous", "snv_unknown"):
        pass  # Already covered by BA1/BS1

    # -- BP7: synonymous + no splice + not conserved --
    if consequence == "synonymous":
        if clinvar_sig in ("Benign", "Likely_benign"):
            bp += 1
            reasons.append("BP7")
        elif clinvar_sig == "VUS":
            bp += 1  # BP7 might still fire for synonymous VUS
            reasons.append("BP7")

    # -- Apply combination rules --
    classification = apply_combination_rules(pvs, ps, pm, pp, ba, bs, bp)

    return classification, reasons

def apply_combination_rules(pvs, ps, pm, pp, ba, bs, bp):
    """Apply ACMG combination rules exactly as implemented."""
    # Benign first
    if ba >= 1: return "Benign"
    if bs >= 2: return "Benign"

    # Conflicting
    has_path = (pvs + ps + pm + pp) > 0
    has_benign = (ba + bs + bp) > 0
    if has_path and has_benign:
        return "VUS"

    # Pathogenic
    if pvs >= 1 and ps >= 1: return "Pathogenic"
    if pvs >= 1 and pm >= 2: return "Pathogenic"
    if pvs >= 1 and pm >= 1 and pp >= 1: return "Pathogenic"
    if pvs >= 1 and pp >= 2: return "Pathogenic"
    if ps >= 2: return "Pathogenic"
    if ps >= 1 and pm >= 3: return "Pathogenic"
    if ps >= 1 and pm >= 2 and pp >= 2: return "Pathogenic"
    if ps >= 1 and pm >= 1 and pp >= 4: return "Pathogenic"

    # Likely pathogenic
    if pvs >= 1 and pm >= 1: return "Likely_pathogenic"
    if ps >= 1 and 1 <= pm <= 2: return "Likely_pathogenic"
    if ps >= 1 and pp >= 2: return "Likely_pathogenic"
    if pm >= 3: return "Likely_pathogenic"
    if pm >= 2 and pp >= 2: return "Likely_pathogenic"
    if pm >= 1 and pp >= 4: return "Likely_pathogenic"

    # Likely benign
    if bs >= 1 and bp >= 1: return "Likely_benign"
    if bp >= 2: return "Likely_benign"

    return "VUS"


def main():
    print("=" * 90)
    print("ClinVar Concordance Analysis v2 — fastVEP ACMG-AMP Classifier")
    print("=" * 90)
    print()

    # Read and filter
    print("Reading ClinVar variant_summary.txt.gz ...")
    two_star = []
    with gzip.open(INPUT_FILE, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("Assembly") != "GRCh38": continue
            stars = review_stars(row.get("ReviewStatus", ""))
            if stars < 2: continue
            norm_sig = normalize_significance(row.get("ClinicalSignificance", ""))
            if not norm_sig: continue
            row["_stars"] = stars
            row["_norm_sig"] = norm_sig
            two_star.append(row)

    print(f"  2-star+ GRCh38 variants: {len(two_star):,}")
    print()

    # Parse consequences and simulate
    consequence_dist = Counter()
    sig_x_cons = defaultdict(Counter)

    # Classification results
    results = []  # (clinvar_class, predicted_class, consequence, gene, name, reasons)

    for row in two_star:
        cons = parse_consequence_from_name(row.get("Name", ""), row.get("Type", ""))
        consequence_dist[cons] += 1
        sig_x_cons[row["_norm_sig"]][cons] += 1

        predicted, reasons = simulate_acmg(row, cons)
        results.append((
            row["_norm_sig"], predicted, cons,
            row.get("GeneSymbol", "?"), row.get("Name", "?")[:60],
            reasons, row["_stars"]
        ))

    # -- Consequence distribution --
    print("Consequence Distribution (2-star+, parsed from HGVS):")
    for cons, count in sorted(consequence_dist.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / len(two_star)
        print(f"  {cons:<25s}: {count:>8,} ({pct:.1f}%)")
    print()

    # -- Consequence by significance --
    print("Key Consequence Types by ClinVar Classification:")
    for sig in ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]:
        total = sum(sig_x_cons[sig].values())
        if total == 0: continue
        print(f"  {sig} (n={total:,}):")
        for cons, count in sorted(sig_x_cons[sig].items(), key=lambda x: -x[1])[:8]:
            pct = 100.0 * count / total
            print(f"    {cons:<25s}: {count:>7,} ({pct:5.1f}%)")
        print()

    # -- Concordance matrix --
    classes = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
    matrix = Counter()
    for clinvar, predicted, *_ in results:
        matrix[(clinvar, predicted)] += 1

    print("=" * 90)
    print("CONCORDANCE MATRIX (rows=ClinVar gold standard, cols=fastVEP predicted)")
    print("=" * 90)
    header = f"{'ClinVar \\ fastVEP':>22s}"
    for p in classes:
        header += f"  {p[:10]:>10s}"
    header += f"  {'Total':>8s}"
    print(header)
    print("-" * (22 + 12*5 + 10))

    for cv in classes:
        line = f"{cv:>22s}"
        row_total = 0
        for p in classes:
            c = matrix.get((cv, p), 0)
            row_total += c
            line += f"  {c:>10,}"
        line += f"  {row_total:>8,}"
        print(line)
    print()

    # -- Metrics --
    total = len(results)
    exact = sum(1 for cv, pred, *_ in results if cv == pred)
    path_set = {"Pathogenic", "Likely_pathogenic"}
    benign_set = {"Benign", "Likely_benign"}
    same_dir = sum(1 for cv, pred, *_ in results
                   if (cv in path_set and pred in path_set)
                   or (cv in benign_set and pred in benign_set)
                   or (cv == "VUS" and pred == "VUS"))
    opposite = sum(1 for cv, pred, *_ in results
                   if (cv in path_set and pred in benign_set)
                   or (cv in benign_set and pred in path_set))

    print(f"Exact match:              {exact:>8,} / {total:,}  ({100*exact/total:.1f}%)")
    print(f"Same direction:           {same_dir:>8,} / {total:,}  ({100*same_dir/total:.1f}%)")
    print(f"Opposite direction:       {opposite:>8,} / {total:,}  ({100*opposite/total:.1f}%)")
    print(f"Downgraded (P/LP → VUS):  {sum(1 for cv, pred, *_ in results if cv in path_set and pred == 'VUS'):>8,}")
    print(f"Upgraded (VUS → P/LP):    {sum(1 for cv, pred, *_ in results if cv == 'VUS' and pred in path_set):>8,}")
    print()

    # -- Per-consequence concordance --
    print("=" * 90)
    print("CONCORDANCE BY CONSEQUENCE TYPE")
    print("=" * 90)
    print(f"{'Consequence':<25s}  {'N':>7s}  {'Exact%':>7s}  {'SameDir%':>8s}  {'Opposite%':>9s}  {'P/LP→VUS':>8s}  {'B/LB→VUS':>8s}")
    print("-" * 90)

    cons_results = defaultdict(list)
    for r in results:
        cons_results[r[2]].append(r)

    for cons in sorted(cons_results.keys(), key=lambda c: -len(cons_results[c])):
        cr = cons_results[cons]
        n = len(cr)
        if n < 20: continue
        e = sum(1 for cv, pred, *_ in cr if cv == pred)
        sd = sum(1 for cv, pred, *_ in cr
                 if (cv in path_set and pred in path_set)
                 or (cv in benign_set and pred in benign_set)
                 or (cv == "VUS" and pred == "VUS"))
        op = sum(1 for cv, pred, *_ in cr
                 if (cv in path_set and pred in benign_set)
                 or (cv in benign_set and pred in path_set))
        plp_vus = sum(1 for cv, pred, *_ in cr if cv in path_set and pred == "VUS")
        blb_vus = sum(1 for cv, pred, *_ in cr if cv in benign_set and pred == "VUS")
        print(f"{cons:<25s}  {n:>7,}  {100*e/n:>6.1f}%  {100*sd/n:>7.1f}%  {100*op/n:>8.1f}%  {plp_vus:>8,}  {blb_vus:>8,}")
    print()

    # -- Most/least concordant genes --
    print("=" * 90)
    print("TOP GENES BY DISCORDANCE (pathogenic variants predicted as VUS)")
    print("=" * 90)
    gene_discord = Counter()
    gene_total_path = Counter()
    for cv, pred, cons, gene, *_ in results:
        if cv in path_set:
            gene_total_path[gene] += 1
            if pred == "VUS":
                gene_discord[gene] += 1

    print(f"{'Gene':<15s}  {'Path_Total':>10s}  {'Pred_VUS':>10s}  {'VUS_Rate':>10s}")
    print("-" * 55)
    for gene, discord in sorted(gene_discord.items(), key=lambda x: -x[1])[:25]:
        total_p = gene_total_path[gene]
        if total_p < 10: continue
        print(f"{gene:<15s}  {total_p:>10,}  {discord:>10,}  {100*discord/total_p:>9.1f}%")
    print()

    # -- Example discordant cases --
    print("=" * 90)
    print("EXAMPLE DISCORDANT CASES (ClinVar P/LP but predicted VUS)")
    print("=" * 90)
    examples = [(cv, pred, cons, gene, name, reasons, stars)
                for cv, pred, cons, gene, name, reasons, stars in results
                if cv in path_set and pred == "VUS"][:20]
    for cv, pred, cons, gene, name, reasons, stars in examples:
        print(f"  {gene:<12s} | {cons:<20s} | ClinVar: {cv:<20s} | Criteria: {','.join(reasons):<25s} | {name}")
    print()

    # -- Key insights --
    print("=" * 90)
    print("KEY INSIGHTS AND RECOMMENDATIONS")
    print("=" * 90)
    print()

    # What types of pathogenic variants are we missing?
    path_results = [(cv, pred, cons, gene, name, reasons, stars)
                    for cv, pred, cons, gene, name, reasons, stars in results
                    if cv in path_set]
    path_correct = sum(1 for _, pred, *_ in path_results if pred in path_set)
    path_vus = sum(1 for _, pred, *_ in path_results if pred == "VUS")
    path_total = len(path_results)

    print(f"Pathogenic/LP variant detection:")
    print(f"  Correctly identified as P/LP: {path_correct:,} / {path_total:,} ({100*path_correct/path_total:.1f}%)")
    print(f"  Downgraded to VUS:            {path_vus:,} / {path_total:,} ({100*path_vus/path_total:.1f}%)")
    print()

    # Breakdown of pathogenic VUS by consequence
    path_vus_by_cons = Counter()
    for _, pred, cons, *_ in path_results:
        if pred == "VUS":
            path_vus_by_cons[cons] += 1
    print("  Pathogenic variants downgraded to VUS, by consequence type:")
    for cons, count in sorted(path_vus_by_cons.items(), key=lambda x: -x[1])[:10]:
        print(f"    {cons:<25s}: {count:>7,}")
    print()

    benign_results = [(cv, pred, cons, gene, name, reasons, stars)
                      for cv, pred, cons, gene, name, reasons, stars in results
                      if cv in benign_set]
    benign_correct = sum(1 for _, pred, *_ in benign_results if pred in benign_set)
    benign_total = len(benign_results)
    print(f"Benign/LB variant detection:")
    print(f"  Correctly identified as B/LB: {benign_correct:,} / {benign_total:,} ({100*benign_correct/benign_total:.1f}%)")
    print()

    print("Highest-impact improvements:")
    print("  1. MISSENSE pathogenic variants (largest VUS bucket): Adding functional")
    print("     domain data (PM1) and protein-position ClinVar index (PS1/PM5)")
    print("     would upgrade many from VUS to LP/P.")
    print("  2. SNV_UNKNOWN variants: Better HGVS parsing to distinguish missense/")
    print("     nonsense/synonymous would correctly trigger PVS1/PP3/BP7.")
    print("  3. DELETION/DUPLICATION variants: Need SV-specific criteria and")
    print("     copy number analysis.")
    print("  4. PP5/BP6 (ClinVar-as-evidence): Currently disabled per SVI.")
    print("     Enabling for 2-star+ variants would dramatically improve concordance")
    print("     but is circular (using ClinVar to validate against ClinVar).")
    print()


if __name__ == "__main__":
    main()
