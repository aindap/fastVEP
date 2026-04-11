#!/usr/bin/env python3
"""Generate manuscript figures for fastVEP."""

import csv
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not installed.")

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.dirname(__file__)

C = {
    'fast': '#6c7aee', 'vep': '#f5426c', 'success': '#10b981',
    'high': '#f5426c', 'moderate': '#f59e0b', 'low': '#3b82f6', 'modifier': '#6b7280',
    'cli': '#6c7aee', 'mid': '#818cf8', 'data': '#34d399', 'core': '#f59e0b', 'sa': '#f472b6',
}
OC = {'yeast': '#10b981', 'drosophila': '#f59e0b', 'arabidopsis': '#8b5cf6', 'mouse': '#ef4444', 'human': '#3b82f6'}
OL = {'yeast': 'Yeast', 'drosophila': 'Drosophila', 'arabidopsis': 'Arabidopsis', 'mouse': 'Mouse', 'human': 'Human'}

def read_csv(f):
    with open(os.path.join(DATA_DIR, f)) as fh:
        return list(csv.DictReader(fh))


def fig1_architecture():
    """Architecture diagram."""
    if not HAS_MPL: return
    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14); ax.set_ylim(0, 9); ax.axis('off')
    def box(x, y, w, h, label, sub, color, fs=11):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.15",facecolor=color,edgecolor='#374151',lw=1.5,alpha=0.9))
        ax.text(x+w/2,y+h/2+0.15,label,ha='center',va='center',fontsize=fs,fontweight='bold',color='white')
        if sub: ax.text(x+w/2,y+h/2-0.2,sub,ha='center',va='center',fontsize=8,color='white',alpha=0.9,style='italic')
    def arr(x1,y1,x2,y2): ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='->',color='#6b7280',lw=1.5))
    ax.text(7,8.6,'fastVEP Architecture',ha='center',fontsize=16,fontweight='bold',color='#1f2937')
    ax.text(7,8.25,'10-crate Cargo workspace',ha='center',fontsize=10,color='#6b7280')
    box(1.5,7,5,0.9,'fastvep-cli','Pipeline, cache, SA builder',C['cli'],13)
    box(7.5,7,5,0.9,'fastvep-web','Production server (axum)',C['cli'],13)
    box(0.3,5.3,2.6,0.9,'fastvep-io','VCF/CSQ/JSON I/O',C['mid'])
    box(3.2,5.3,2.6,0.9,'fastvep-hgvs','HGVSg/c/p',C['mid'])
    box(6.1,5.3,2.8,0.9,'fastvep-consequence','SNV/indel/SV',C['mid'])
    box(9.2,5.3,2.2,0.9,'fastvep-filter','Filter engine',C['mid'])
    box(11.7,5.3,2,0.9,'fastvep-sa','Annotations',C['sa'])
    box(0.5,3.5,3.5,0.9,'fastvep-genome','Transcript, Exon, Gene',C['data'])
    box(4.5,3.5,5,0.9,'fastvep-cache','GFF3, FASTA mmap, cache',C['data'])
    box(10,3.5,3.5,0.9,'fastvep-sa (fastSA)','ClinVar, gnomAD, ...',C['sa'])
    box(3,1.7,8,0.9,'fastvep-core','Consequence (49 SO), Allele, Impact, VariantType',C['core'],12)
    for x in [1.6,4.5,7.5,10.3,12.7]: arr(4,7,x,6.2); arr(10,7,x,6.2)
    for x in [2.25,7,11.75]: arr(x,5.3,x,4.4); arr(x,3.5,7,2.6)
    ax.text(7,0.8,'17,966 LOC  |  175 tests  |  3.3 MB binary  |  fastVEP.org',ha='center',fontsize=10,color='#6b7280',bbox=dict(boxstyle='round,pad=0.3',facecolor='#f3f4f6',edgecolor='#d1d5db'))
    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig1_architecture.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig1_architecture"); plt.close()


def fig2_comparison():
    """THE figure: fastVEP vs VEP with all data on one plot."""
    h2h = read_csv('vep_comparison.csv')
    orgs = read_csv('organism_comparison.csv')
    if not HAS_MPL: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

    # === Panel A: Wall-clock time (log-log) — both tools on same axes ===
    # VEP h2h points
    hv = [int(d['variants']) for d in h2h]
    ft = [float(d['fastvep_sec']) for d in h2h]
    vt = [float(d['vep_sec']) for d in h2h]
    ax1.plot(hv, ft, 'o-', color=C['fast'], lw=2.5, ms=10, label='fastVEP (chr22)', zorder=3)
    ax1.plot(hv, vt, 's-', color=C['vep'], lw=2.5, ms=10, label='Ensembl VEP v115', zorder=3)

    # VEP extrapolation (dashed, VEP can't do more)
    ax1.plot([50000, 4048342], [206.13, 206.13*(4048342/50000)], 's--', color=C['vep'], lw=1.5, ms=6, alpha=0.3)
    ax1.annotate('VEP: cannot complete\nfull genome (~4.6 hrs est.)', xy=(4048342, 206.13*(4048342/50000)),
                 xytext=(-30, -25), textcoords='offset points', fontsize=8, color=C['vep'], style='italic', ha='right',
                 arrowprops=dict(arrowstyle='->', color=C['vep'], lw=1, ls='--'))

    # fastVEP multi-organism points (all on same plot!)
    for d in orgs:
        o = d['organism']
        v, t = int(d['variants']), float(d['time_sec'])
        ax1.scatter(v, t, c=OC.get(o,'#999'), s=180, marker='D', alpha=0.9, edgecolors='white', lw=1.5, zorder=4)
        label = f"{OL.get(o,o)}\n{v/1e6:.1f}M" if v > 500000 else f"{OL.get(o,o)}\n{v/1e3:.0f}K"
        ax1.annotate(label, xy=(v,t), xytext=(10, 6), textcoords='offset points', fontsize=7.5, fontweight='bold', color='#333')

    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_xlabel('Variants', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. fastVEP vs Ensembl VEP: Annotation Time', fontsize=13, fontweight='bold')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,_: f'{x/1e6:.0f}M' if x>=1e6 else f'{x/1e3:.0f}K'))
    ax1.grid(True, alpha=0.15)

    # Custom legend
    legend_h = [
        plt.Line2D([0],[0],marker='o',color=C['fast'],lw=2,ms=8,label='fastVEP (chr22 h2h)'),
        plt.Line2D([0],[0],marker='s',color=C['vep'],lw=2,ms=8,label='Ensembl VEP v115'),
        plt.Line2D([0],[0],marker='D',color='#999',lw=0,ms=8,label='fastVEP (full organisms)'),
    ]
    ax1.legend(handles=legend_h, fontsize=9, loc='upper left')

    # === Panel B: Speedup + throughput ===
    # Speedup bars (left y-axis)
    speedups = [float(d['speedup']) for d in h2h]
    labels = [f"{int(d['variants'])/1000:.0f}K" for d in h2h]
    x = range(len(labels))
    bars = ax2.bar(x, speedups, color=C['fast'], alpha=0.85, width=0.6)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=11)
    ax2.set_xlabel('Variants (GIAB HG002 chr22)', fontsize=12)
    ax2.set_ylabel('Speedup over Ensembl VEP', fontsize=12, color=C['fast'])
    ax2.set_title('B. Speedup: fastVEP / VEP', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.15, axis='y')

    for bar, sp in zip(bars, speedups):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                 f'{sp:.0f}x', ha='center', fontsize=12, fontweight='bold', color=C['fast'])

    # VEP throughput degradation line (secondary y-axis)
    ax2b = ax2.twinx()
    vep_vps = [int(d['vep_vps']) for d in h2h]
    ax2b.plot(x, vep_vps, 's--', color=C['vep'], lw=2, ms=8, alpha=0.7, label='VEP throughput')
    ax2b.set_ylabel('VEP throughput (v/s)', fontsize=11, color=C['vep'])
    ax2b.tick_params(axis='y', labelcolor=C['vep'])

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig2_comparison.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig2_comparison"); plt.close()


def fig3_concordance():
    """100% VEP concordance."""
    data = read_csv('vep_concordance.csv')
    if not HAS_MPL: return
    fields = [d['field'] for d in data]
    accuracy = [float(d['accuracy']) for d in data]
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = [C['success'] if a==100.0 else C['vep'] for a in accuracy]
    bars = ax.barh(range(len(fields)), accuracy, color=colors, alpha=0.85)
    ax.set_yticks(range(len(fields)))
    ax.set_yticklabels(fields, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Concordance with Ensembl VEP v115.1 (%)', fontsize=12)
    ax.set_title('Field-Level Accuracy: fastVEP vs Ensembl VEP\n(2,340 shared transcript-allele pairs, 173 variants)', fontsize=13, fontweight='bold')
    ax.set_xlim(95, 101); ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=100, color=C['success'], linestyle='--', alpha=0.5, lw=1)
    for bar, acc in zip(bars, accuracy):
        ax.text(bar.get_width()-0.3, bar.get_y()+bar.get_height()/2, f'{acc:.1f}%', va='center', ha='right', fontsize=9, color='white', fontweight='bold')
    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig3_concordance.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig3_concordance"); plt.close()


def fig4_consequences():
    """Consequence distribution."""
    data = read_csv('consequence_distribution.csv')
    if not HAS_MPL: return
    consequences = [d['consequence'] for d in data]
    counts = [int(d['count']) for d in data]
    impact_map = {'splice_acceptor_variant':'HIGH','splice_donor_variant':'HIGH','stop_gained':'HIGH','frameshift_variant':'HIGH',
                  'missense_variant':'MODERATE','inframe_insertion':'MODERATE','inframe_deletion':'MODERATE',
                  'splice_region_variant':'LOW','synonymous_variant':'LOW','splice_polypyrimidine_tract_variant':'LOW',
                  'splice_donor_5th_base_variant':'LOW','splice_donor_region_variant':'LOW'}
    colors = [C.get(impact_map.get(c,'MODIFIER').lower(), C['modifier']) for c in consequences]
    fig, ax = plt.subplots(figsize=(11, 7))
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences))); ax.set_yticklabels(consequences, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis(); ax.set_xlabel('Annotation count', fontsize=12)
    ax.set_title('Predicted Consequence Distribution\n(GIAB HG002 full WGS: 4.05M variants, 50.1M annotations)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x'); ax.set_xscale('log')
    for bar, count in zip(bars, counts):
        ax.text(bar.get_width()*1.1, bar.get_y()+bar.get_height()/2, f'{count:,}', va='center', fontsize=9, color='#555')
    legend_el = [mpatches.Patch(facecolor=C[k],alpha=0.85,label=l) for k,l in [('high','HIGH'),('moderate','MODERATE'),('low','LOW'),('modifier','MODIFIER')]]
    ax.legend(handles=legend_el, title='Impact', loc='lower right', fontsize=9)
    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig4_consequences.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig4_consequences"); plt.close()


def fig5_organisms():
    """Multi-organism throughput summary bar chart."""
    data = read_csv('organism_comparison.csv')
    if not HAS_MPL: return

    orgs = [d['organism'] for d in data]
    vps = [int(d['variants_per_sec']) for d in data]
    variants = [int(d['variants']) for d in data]
    times = [float(d['time_sec']) for d in data]
    transcripts = [int(d['transcripts']) for d in data]
    colors = [OC.get(o,'#999') for o in orgs]

    # Sort by throughput
    idx = sorted(range(len(orgs)), key=lambda i: vps[i])
    labels = [f"{OL.get(orgs[i],orgs[i])}\n({transcripts[i]:,} tr., {variants[i]/1e6:.1f}M var.)" for i in idx]

    fig, ax = plt.subplots(figsize=(12, 5.5))
    bars = ax.barh(range(len(idx)), [vps[i] for i in idx], color=[colors[i] for i in idx], alpha=0.85, height=0.6)
    ax.set_yticks(range(len(idx))); ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Throughput (variants/sec)', fontsize=12)
    ax.set_title('Gold-Standard Annotation Throughput by Organism\n(complete VCFs, full Ensembl genome annotations, HGVS enabled)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.15, axis='x')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,_: f'{x/1000:.0f}K'))

    for bar, i in zip(bars, idx):
        ax.text(bar.get_width()+1000, bar.get_y()+bar.get_height()/2,
                f'{vps[i]:,} v/s  ({times[i]:.0f}s)', va='center', fontsize=10, fontweight='bold', color='#333')

    plt.tight_layout()
    for ext in ['png','pdf']: plt.savefig(os.path.join(FIG_DIR, f'fig5_organisms.{ext}'), dpi=300, bbox_inches='tight')
    print("Saved fig5_organisms"); plt.close()


if __name__ == '__main__':
    print("Generating fastVEP manuscript figures...\n")
    fig1_architecture()
    fig2_comparison()
    fig3_concordance()
    fig4_consequences()
    fig5_organisms()
    print("\nDone.")
