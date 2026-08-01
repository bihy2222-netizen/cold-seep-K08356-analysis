from pathlib import Path
import sys, re, csv

hmm_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
input48 = Path(sys.argv[3])
allmag = Path(sys.argv[4])

hmmres = out_dir / '01_hmm_result'
tables = out_dir / '03_tables'
newdir = out_dir / '04_new_candidates'
tables.mkdir(parents=True, exist_ok=True); newdir.mkdir(parents=True, exist_ok=True)

def fasta_ids(path):
    ids=[]
    for line in path.open(errors='ignore'):
        if line.startswith('>'):
            ids.append(line[1:].split()[0])
    return ids

def norm_id(x):
    # DRAM genes.faa uses headers like BIN_BIN-k141_contig_gene; tree tips keep only BIN-k141_contig_gene.
    prefixes = ('C1','C2','C3','R2111','S1','S13','S14','S15','SQ','SY')
    candidates=[]
    for i in range(len(x)):
        if i > 0 and x[i-1] not in {'_', '|', ' '} :
            continue
        if not x.startswith(prefixes, i):
            continue
        sub=x[i:]
        m=re.match(r'((?:C\d|R2111|S\d+|SQ|SY)[^\s|]*?_bin\d+-k141_\d+_\d+)$', sub)
        if m:
            candidates.append(m.group(1))
    if candidates:
        return candidates[-1]
    return x

def read_tbl(path):
    best={}
    if not path.exists(): return best
    for line in path.open(errors='ignore'):
        if not line.strip() or line.startswith('#'): continue
        p=line.split()
        if len(p) < 6: continue
        target=p[0]
        try:
            e=float(p[4]); score=float(p[5]); bias=float(p[6]) if len(p)>6 else 0.0
        except ValueError:
            continue
        old=best.get(target)
        if old is None or score > old['bitscore']:
            best[target]={'evalue':e,'bitscore':score,'bias':bias}
    return best

combined=read_tbl(hmmres/'combined_48.tbl')
idr=read_tbl(hmmres/'idrA_48.tbl')
aio=read_tbl(hmmres/'aioA_48.tbl')
ids48=fasta_ids(input48)
with (tables/'HMM_scores_48.tsv').open('w', newline='') as fo:
    w=csv.writer(fo, delimiter='\t')
    w.writerow(['protein_id','combined_bitscore','combined_evalue','idrA_bitscore','idrA_evalue','aioA_bitscore','aioA_evalue','idrA_minus_aioA','hmm_score_note'])
    for pid in ids48:
        c=combined.get(pid,{}); i=idr.get(pid,{}); a=aio.get(pid,{})
        ib=i.get('bitscore'); ab=a.get('bitscore')
        diff='' if ib is None or ab is None else round(ib-ab,3)
        if ib is None and ab is None: note='no_idrA_or_aioA_hit'
        elif ib is not None and (ab is None or ib>ab): note='idrA_model_higher'
        elif ab is not None and (ib is None or ab>ib): note='aioA_model_higher'
        else: note='idrA_aioA_tie_or_unclear'
        w.writerow([pid,c.get('bitscore',''),c.get('evalue',''),i.get('bitscore',''),i.get('evalue',''),a.get('bitscore',''),a.get('evalue',''),diff,note])

allhits=read_tbl(hmmres/'all_MAG_combined.tbl')
with (tables/'all_MAG_T640_hits.tsv').open('w', newline='') as fo:
    w=csv.writer(fo, delimiter='\t')
    w.writerow(['protein_id','normalized_id','combined_bitscore','combined_evalue'])
    for pid,row in sorted(allhits.items(), key=lambda kv:(-kv[1]['bitscore'], kv[0])):
        w.writerow([pid,norm_id(pid),row['bitscore'],row['evalue']])

known_norm={norm_id(x) for x in ids48}
new_ids=[pid for pid in allhits if norm_id(pid) not in known_norm]
(tables/'new_T640_hit_ids.txt').write_text('\n'.join(sorted(new_ids)) + ('\n' if new_ids else ''))

# extract new candidates from all MAG FASTA
want=set(new_ids); keep=False; nnew=0
with allmag.open(errors='ignore') as fi, (newdir/'new_T640_hits.faa').open('w') as fo:
    for line in fi:
        if line.startswith('>'):
            pid=line[1:].split()[0]
            keep=pid in want
            if keep: nnew += 1
        if keep: fo.write(line)

with (tables/'final_evidence_template.tsv').open('w', newline='') as fo:
    w=csv.writer(fo, delimiter='\t')
    w.writerow(['protein_id','normalized_id','source','combined_bitscore','idrA_bitscore','aioA_bitscore','idrA_minus_aioA','tree_clade','gene_neighborhood','final_class','confidence','notes'])
    for pid in ids48:
        i=idr.get(pid,{}); a=aio.get(pid,{}); c=combined.get(pid,{})
        ib=i.get('bitscore'); ab=a.get('bitscore')
        diff='' if ib is None or ab is None else round(ib-ab,3)
        w.writerow([pid,norm_id(pid),'original_48',c.get('bitscore',''),i.get('bitscore',''),a.get('bitscore',''),diff,'','','','',''])
    for pid in sorted(new_ids):
        c=allhits[pid]
        w.writerow([pid,norm_id(pid),'new_all_MAG_T640_hit',c.get('bitscore',''),'','','','','','','',''])

summary = []
summary.append(f'input_48_candidates\t{len(ids48)}')
summary.append(f'combined_48_hits\t{len(combined)}')
summary.append(f'idrA_48_hits\t{len(idr)}')
summary.append(f'aioA_48_hits\t{len(aio)}')
summary.append(f'all_MAG_T640_hits\t{len(allhits)}')
summary.append(f'new_T640_hits_not_in_48_by_normalized_id\t{len(new_ids)}')
(tables/'run_summary.txt').write_text('\n'.join(summary)+'\n')
print('\n'.join(summary))
