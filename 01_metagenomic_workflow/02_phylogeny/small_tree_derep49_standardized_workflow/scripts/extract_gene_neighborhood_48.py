from pathlib import Path
import re, csv
from collections import defaultdict

cand_faa = Path('input/48_candidates.faa')
genes_faa = Path('/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/all_50_10_dRep/dereplicated_genomes/bin_DRAM/annotation/genes.faa')
out_long = Path('result/03_tables/gene_neighborhood_48_long.tsv')
out_summary = Path('result/03_tables/gene_neighborhood_48_summary.tsv')

def parse_target_id(pid):
    m = re.match(r'(.+?_bin\d+)-k141_(\d+)_(\d+)$', pid)
    if not m:
        return None
    bin_id, contig_num, gene_num = m.groups()
    contig_core = f'{bin_id}-k141_{contig_num}'
    return bin_id, contig_core, int(gene_num)

def parse_gene_header(line):
    raw = line[1:].rstrip('\n')
    first = raw.split()[0]
    m = re.match(r'(.+?_bin\d+)_(.+?_bin\d+-k141_(\d+)_(\d+))$', first)
    if not m:
        return None
    bin_prefix, full_core, contig_num, gene_num = m.groups()
    contig_core = re.sub(r'^(.+?_bin\d+)_', '', full_core)
    # For headers with repeated bin prefix, contig_core should be BIN-k141_n.
    m2 = re.match(r'(.+?_bin\d+-k141_\d+)_\d+$', contig_core)
    if not m2:
        return None
    contig_id = m2.group(1)
    annot = raw.split(' ', 1)[1] if ' ' in raw else ''
    return {
        'raw_id': first,
        'normalized_id': contig_core,
        'bin_id': bin_prefix,
        'contig_id': contig_id,
        'gene_index': int(gene_num),
        'annotation': annot,
    }

def has_term(text, terms):
    low = text.lower()
    return any(t.lower() in low for t in terms)

targets=[]
for line in cand_faa.open():
    if line.startswith('>'):
        pid=line[1:].split()[0]
        parsed=parse_target_id(pid)
        if parsed:
            targets.append({'target_id': pid, 'bin_id': parsed[0], 'contig_id': parsed[1], 'target_gene_index': parsed[2]})

target_contigs={t['contig_id'] for t in targets}
by_contig=defaultdict(list)
for line in genes_faa.open(errors='ignore'):
    if line.startswith('>'):
        rec=parse_gene_header(line)
        if rec and rec['contig_id'] in target_contigs:
            by_contig[rec['contig_id']].append(rec)
for rows in by_contig.values():
    rows.sort(key=lambda r: r['gene_index'])

long_rows=[]; summaries=[]
for t in targets:
    rows=by_contig.get(t['contig_id'], [])
    window=[r for r in rows if abs(r['gene_index']-t['target_gene_index']) <= 10]
    target_row=next((r for r in rows if r['gene_index']==t['target_gene_index']), None)
    flags={
        'idrB_present': False,
        'idrP1_present': False,
        'idrP2_present': False,
        'aioB_present': False,
        'aio_operon_like': False,
    }
    neighbor_terms=[]
    for r in window:
        ann = r['annotation']
        if r['gene_index'] == t['target_gene_index']:
            role='target'
        elif r['gene_index'] < t['target_gene_index']:
            role='upstream'
        else:
            role='downstream'
        rel=r['gene_index']-t['target_gene_index']
        gene_label=''
        if has_term(ann, ['idrB']): flags['idrB_present']=True; gene_label='idrB'
        if has_term(ann, ['idrP1','idrP 1','idrP_1']): flags['idrP1_present']=True; gene_label='idrP1'
        if has_term(ann, ['idrP2','idrP 2','idrP_2']): flags['idrP2_present']=True; gene_label='idrP2'
        if has_term(ann, ['aioB','arsenite oxidase small subunit']): flags['aioB_present']=True; gene_label='aioB'
        if has_term(ann, ['aioA','aioB','arsenite oxidase','molybdopterin oxidoreductase','4fe-4s','tat']):
            neighbor_terms.append(f'{rel}:{ann[:90]}')
        long_rows.append({
            'target_id': t['target_id'], 'bin_id': t['bin_id'], 'contig_id': t['contig_id'],
            'target_gene_index': t['target_gene_index'], 'neighbor_gene_index': r['gene_index'],
            'relative_position': rel, 'role': role, 'neighbor_raw_id': r['raw_id'],
            'neighbor_normalized_id': r['normalized_id'], 'putative_gene_label': gene_label,
            'annotation': ann,
        })
    flags['aio_operon_like'] = flags['aioB_present']
    if flags['idrB_present'] and flags['idrP1_present'] and flags['idrP2_present']:
        neigh = 'complete idrB-idrP1-idrP2-like neighborhood'
    elif flags['idrB_present'] or flags['idrP1_present'] or flags['idrP2_present']:
        neigh = 'partial Idr-associated neighborhood'
    elif flags['aioB_present']:
        neigh = 'aioB/aio-like neighborhood'
    elif window:
        neigh = 'no idrB/idrP1/idrP2/aioB detected in DRAM +/-10 annotations'
    else:
        neigh = 'target contig not found in DRAM genes.faa'
    summaries.append({
        'target_id': t['target_id'], 'bin_id': t['bin_id'], 'contig_id': t['contig_id'],
        'target_gene_index': t['target_gene_index'], 'n_genes_in_window': len(window),
        'target_annotation': target_row['annotation'] if target_row else '',
        **{k: ('Yes' if v else 'No') for k,v in flags.items()},
        'gene_neighborhood': neigh,
        'neighborhood_notes': '; '.join(neighbor_terms[:8]),
    })

with out_long.open('w', newline='') as f:
    fields=list(long_rows[0].keys()) if long_rows else ['target_id']
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(long_rows)
with out_summary.open('w', newline='') as f:
    fields=list(summaries[0].keys()) if summaries else ['target_id']
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(summaries)
print('targets', len(targets))
print('summary', out_summary)
print('long', out_long)
