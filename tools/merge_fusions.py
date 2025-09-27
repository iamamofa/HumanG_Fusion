#!/usr/bin/env python3
"""
merge_fusions.py - very small merger for STAR-Fusion and Arriba (and fusioncatcher)
Usage:
  merge_fusions.py --starfusion <file> --arriba <file> --fusioncatcher <glob> --gtf <gtf> --out <out>
"""
import argparse, csv, gzip, glob
from collections import defaultdict

def read_starfusion(path):
    calls = []
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for r in reader:
                calls.append({
                    'caller':'STAR-Fusion',
                    'fusion_name': r.get('#FusionName') or r.get('FusionName') or r.get('LeftGene')+'--'+r.get('RightGene'),
                    'left_gene': r.get('LeftGene'),
                    'right_gene': r.get('RightGene'),
                    'split_reads': r.get('JunctionReads') or r.get('SplitReads'),
                    'spanning_reads': r.get('SpanningFragCount') or r.get('SpanningReads')
                })
    except Exception:
        pass
    return calls

def read_arriba(path):
    calls=[]
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'): continue
                cols=line.strip().split('\t')
                # arriba columns: chr1, pos1, strand1, gene1, chr2, pos2, strand2, gene2, ...
                gene1=cols[3]
                gene2=cols[7]
                support = cols[8] if len(cols)>8 else ''
                calls.append({
                    'caller':'Arriba',
                    'fusion_name': gene1+'--'+gene2,
                    'left_gene': gene1,
                    'right_gene': gene2,
                    'support': support
                })
    except Exception:
        pass
    return calls

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--starfusion', default=None)
    p.add_argument('--arriba', default=None)
    p.add_argument('--fusioncatcher', default=None)
    p.add_argument('--gtf', required=True)
    p.add_argument('--out', required=True)
    args=p.parse_args()

    merged=defaultdict(lambda: {'callers':set(), 'left':None, 'right':None, 'support':[]})
    if args.starfusion:
        for c in read_starfusion(args.starfusion):
            key=(c['left_gene'], c['right_gene'])
            merged[key]['callers'].add('STAR-Fusion'); merged[key]['left']=c['left_gene']; merged[key]['right']=c['right_gene']; merged[key]['support'].append(c.get('split_reads') or '')
    if args.arriba:
        for c in read_arriba(args.arriba):
            key=(c['left_gene'], c['right_gene'])
            merged[key]['callers'].add('Arriba'); merged[key]['left']=c['left_gene']; merged[key]['right']=c['right_gene']; merged[key]['support'].append(c.get('support') or '')

    # FusionCatcher parsing left as exercise - here we skip for brevity

    # write TSV
    with open(args.out,'w') as outfh:
        outfh.write("left_gene\tright_gene\tcallers\tsupport_count\tsupports\n")
        for k,v in merged.items():
            outfh.write(f"{k[0]}\t{k[1]}\t{';'.join(sorted(v['callers']))}\t{len(v['support'])}\t{';'.join(v['support'])}\n")

if __name__=='__main__':
    main()
