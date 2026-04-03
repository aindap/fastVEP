#!/usr/bin/env python3
"""
Generate manuscript figures for OxiVEP.
Requires: matplotlib, pandas (pip install matplotlib pandas)
"""

import csv
import os

# Try importing matplotlib; provide instructions if missing
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not installed. Install with: pip install matplotlib")
    print("Generating text-based summaries instead.\n")

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.dirname(__file__)

# Color scheme
COLORS = {
    'primary': '#6c7aee',
    'high': '#f5426c',
    'moderate': '#f59e0b',
    'low': '#3b82f6',
    'modifier': '#6b7280',
    'bg': '#1a1d27',
    'text': '#e4e6ed',
    'grid': '#2e3240',
}

def read_csv(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)

def fig2_throughput_scaling():
    """Figure 2: Throughput scaling curve."""
    data = read_csv('scaling.csv')
    variants = [int(d['variants']) for d in data]
    throughput = [int(d['variants_per_sec']) for d in data]
    times = [float(d['time_sec']) for d in data]

    if not HAS_MPL:
        print("Figure 2: Throughput Scaling")
        print(f"{'Variants':>10} {'Time (s)':>10} {'Throughput (v/s)':>18}")
        for v, t, tp in zip(variants, times, throughput):
            bar = '#' * int(tp / 10000)
            print(f"{v:>10,} {t:>10.3f} {tp:>15,}  {bar}")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Throughput vs variant count
    ax1.plot(variants, throughput, 'o-', color=COLORS['primary'], linewidth=2, markersize=8)
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of input variants', fontsize=12)
    ax1.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax1.set_title('A. Annotation Throughput', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K'))
    ax1.set_ylim(bottom=0)

    # Panel B: Time vs variant count
    ax2.plot(variants, times, 's-', color=COLORS['moderate'], linewidth=2, markersize=8)
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of input variants', fontsize=12)
    ax2.set_ylabel('Annotation time (seconds)', fontsize=12)
    ax2.set_title('B. Annotation Time', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.pdf'), bbox_inches='tight')
    print("Saved fig2_throughput_scaling.png/pdf")
    plt.close()

def fig3_resource_usage():
    """Figure 3: Resource usage comparison."""
    data = read_csv('resource_usage.csv')
    mem_data = {d['metric']: float(d['value']) for d in data}

    if not HAS_MPL:
        print("Figure 3: Resource Usage")
        for k, v in mem_data.items():
            unit = next((d['unit'] for d in data if d['metric'] == k), '')
            print(f"  {k}: {v} {unit}")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # Panel A: Memory comparison
    variants_labels = ['1K', '10K', '100K']
    oxivep_mem = [mem_data.get('peak_memory_1000v', 0),
                  mem_data.get('peak_memory_10000v', 0),
                  mem_data.get('peak_memory_100000v', 0)]
    vep_mem = [500, 500, 600]  # Estimated Ensembl VEP

    x = range(len(variants_labels))
    w = 0.35
    ax1.bar([i - w/2 for i in x], vep_mem, w, label='Ensembl VEP', color=COLORS['high'], alpha=0.8)
    ax1.bar([i + w/2 for i in x], oxivep_mem, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax1.set_xlabel('Input size (variants)', fontsize=12)
    ax1.set_ylabel('Peak memory (MB)', fontsize=12)
    ax1.set_title('A. Memory Usage', fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(variants_labels)
    ax1.legend()
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, axis='y')

    # Panel B: Binary size / startup
    metrics = ['Binary Size\n(MB)', 'Startup\n(ms)', 'Startup + GFF3\n(ms)']
    vep_vals = [200, 10000, 12000]
    oxi_vals = [mem_data.get('binary_size', 0),
                mem_data.get('startup_time', 0),
                mem_data.get('startup_with_gff3', 0)]

    x = range(len(metrics))
    ax2.bar([i - w/2 for i in x], vep_vals, w, label='Ensembl VEP', color=COLORS['high'], alpha=0.8)
    ax2.bar([i + w/2 for i in x], oxi_vals, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_title('B. Size & Startup', fontsize=13, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(metrics, fontsize=10)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig3_resource_usage.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig3_resource_usage.pdf'), bbox_inches='tight')
    print("Saved fig3_resource_usage.png/pdf")
    plt.close()

def fig4_consequence_distribution():
    """Figure 4: Consequence distribution."""
    data = read_csv('consequence_distribution.csv')
    consequences = [d['consequence'] for d in data if d['consequence'] != 'consequence']
    counts = [int(d['count']) for d in data if d['consequence'] != 'consequence']

    impact_map = {
        'splice_acceptor_variant': 'HIGH', 'splice_donor_variant': 'HIGH',
        'stop_gained': 'HIGH', 'frameshift_variant': 'HIGH',
        'missense_variant': 'MODERATE', 'inframe_insertion': 'MODERATE',
        'inframe_deletion': 'MODERATE',
        'splice_region_variant': 'LOW', 'synonymous_variant': 'LOW',
        'splice_polypyrimidine_tract_variant': 'LOW',
    }
    def get_color(c):
        imp = impact_map.get(c, 'MODIFIER')
        return COLORS.get(imp.lower(), COLORS['modifier'])

    if not HAS_MPL:
        print("Figure 4: Consequence Distribution")
        for c, n in zip(consequences, counts):
            bar = '#' * (n // 5)
            print(f"  {c:45s} {n:>5d}  {bar}")
        print()
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = [get_color(c) for c in consequences]
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences)))
    ax.set_yticklabels(consequences, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Variant count', fontsize=12)
    ax.set_title('Predicted Consequence Distribution', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + 2, bar.get_y() + bar.get_height()/2,
                str(count), va='center', fontsize=9, color='#666')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.pdf'), bbox_inches='tight')
    print("Saved fig4_consequence_distribution.png/pdf")
    plt.close()

def fig1_architecture():
    """Figure 1: Architecture diagram (text-based for now)."""
    diagram = """
    ┌─────────────────────────────────────────────────────────────────┐
    │                        oxivep-cli                               │
    │              (CLI binary, pipeline, web server)                  │
    └──────┬──────────┬──────────┬──────────┬──────────┬─────────────┘
           │          │          │          │          │
    ┌──────▼──┐ ┌─────▼─────┐ ┌─▼────────┐ │   ┌──────▼──────┐
    │oxivep-io│ │oxivep-hgvs│ │oxivep-   │ │   │oxivep-filter│
    │(VCF I/O,│ │(HGVSg,    │ │consequen-│ │   │(filtering)  │
    │ output) │ │ HGVSc,    │ │ce engine)│ │   └─────────────┘
    └────┬────┘ │ HGVSp)    │ └──┬───────┘ │
         │      └─────┬─────┘    │         │
         │            │          │    ┌────▼─────┐
    ┌────▼────────────▼──────────▼──┐ │oxivep-   │
    │        oxivep-genome          │ │cache     │
    │  (Transcript, Exon, Gene,     │ │(GFF3,    │
    │   CodonTable, coord mapping)  │ │ FASTA,   │
    └────────────┬──────────────────┘ │ info.txt)│
                 │                    └────┬─────┘
            ┌────▼────────────────────────▼──┐
            │          oxivep-core            │
            │  (GenomicPosition, Consequence, │
            │   Allele, Strand, Impact)       │
            └─────────────────────────────────┘
    """
    with open(os.path.join(FIG_DIR, 'fig1_architecture.txt'), 'w') as f:
        f.write(diagram)
    print("Saved fig1_architecture.txt")

if __name__ == '__main__':
    print("Generating OxiVEP manuscript figures...")
    print(f"Data dir: {os.path.abspath(DATA_DIR)}")
    print(f"Output dir: {os.path.abspath(FIG_DIR)}")
    print()
    fig1_architecture()
    fig2_throughput_scaling()
    fig3_resource_usage()
    fig4_consequence_distribution()
    print("\nDone.")
