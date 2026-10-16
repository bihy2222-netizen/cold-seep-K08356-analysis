from pathlib import Path
import csv,json,collections,shutil
O=Path(__file__).resolve().parent;P=O/'inputs';T=O/'tables'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def save(name,headers,rows):
 with (T/name).open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(headers);w.writerows(rows)
G=[('K00200','fwdA/fmdA','一碳转化','甲酰甲烷呋喃脱氢酶A亚基'),('K00201','fwdB/fmdB','一碳转化','甲酰甲烷呋喃脱氢酶B亚基'),('K00202','fwdC/fmdC','一碳转化','甲酰甲烷呋喃脱氢酶C亚基'),('K00203','fwdD/fmdD','一碳转化','甲酰甲烷呋喃脱氢酶D亚基'),('K00205','fwdF/fmdF','一碳转化及电子传递','4Fe-4S铁氧还蛋白'),('K00319','mtd','一碳转化','亚甲基四氢甲烷蝶呤脱氢酶'),('K00320','mer','一碳转化','亚甲基四氢甲烷蝶呤还原酶'),('K00399','mcrA','甲烷生成或活化','甲基辅酶M还原酶α亚基'),('K00401','mcrB','甲烷生成或活化','甲基辅酶M还原酶β亚基'),('K00402','mcrG','甲烷生成或活化','甲基辅酶M还原酶γ亚基')]
G += [(f'K{577+i:05d}','mtr'+c,'甲基转移','四氢甲烷蝶呤S-甲基转移酶'+c+'亚基') for i,c in enumerate('ABCDEFGH')]
G += [('K00672','ftr','一碳转化','甲酰甲烷呋喃—四氢甲烷蝶呤甲酰基转移酶'),('K01499','mch','一碳转化','次甲基四氢甲烷蝶呤环水解酶'),('K03388','hdrA2','辅酶再生及电子传递','异二硫键还原酶A2亚基'),('K03389','hdrB2','辅酶再生及电子传递','异二硫键还原酶B2亚基'),('K03390','hdrC2','辅酶再生及电子传递','异二硫键还原酶C2亚基')]
byko={x[0]:x for x in G};ks=list(byko);ko_by_gene={x[1]:x[0] for x in G}
pathgenes=['mcrA','mcrB','mcrG']+['mtr'+c for c in 'ABCDEFGH']+['mer','mtd','mch','ftr']+['fwd'+c+'/fmd'+c for c in 'ABCDF']+['hdrA2','hdrB2','hdrC2']
refgenes=['mcrA','mcrG','mcrB']+['mtr'+c for c in 'HGFABCDE']+['mer','mtd','mch','ftr']+['fwd'+c+'/fmd'+c for c in 'CABDF']+['hdrA2','hdrB2','hdrC2']
order=[ko_by_gene[x] for x in pathgenes];reforder=[ko_by_gene[x] for x in refgenes];assert len(set(order))==len(order)==23
meta=read(T/'cold_seep_all808_taxonomy_quality.tsv');md={x['MAG_ID']:x for x in meta};assert len(md)==808
selected=read(P/'tables/mcrA_positive_dereplicated_MAGs.tsv');sel={x['selected_positive_MAG']:x for x in selected};cls={x['MAG_ID']:x for x in read(P/'tables/selected_MAG_McrA_classification.tsv')}
for m,x in sel.items():
 if m not in md:
  md[m]={'MAG_ID':m,'classification':cls[m]['host_taxonomy'],'completeness':'NA','contamination':'NA','dRep_cluster':x['dRep_cluster']}
counts=collections.defaultdict(set);hits={}
for x in read(P/'tables/representative_KO_hits.tsv'):
 if x['KO'] in byko:hits[x['MAG_ID'],x['protein_ID'],x['KO']]=[x['MAG_ID'],x['protein_ID'],x['KO'],byko[x['KO']][1],'Kofam严格命中',x['score'],x['E_value'],'原808代表']
for x in read(P/'tables/representative_HMM_hits.tsv'):
 if x['KO'] in byko and x['screen_status']=='strict':hits.setdefault((x['MAG_ID'],x['protein_ID'],x['KO']),[x['MAG_ID'],x['protein_ID'],x['KO'],byko[x['KO']][1],'HMM严格命中',x['score'],x['domain_i_Evalue'],'原808代表'])
for x in read(P/'tables/extra_selected_MAG_pathway_HMM_hits.tsv'):
 if x['KO'] in byko:hits[x['MAG_ID'],x['protein_ID'],x['KO']]=[x['MAG_ID'],x['protein_ID'],x['KO'],byko[x['KO']][1],'补选MAG HMM严格命中',x['score'],'NA','补选阳性代表']
for m,x in cls.items():
 key=m,x['McrA_protein_ID'],'K00399';hits.setdefault(key,[m,x['McrA_protein_ID'],'K00399','mcrA',x['McrA_evidence_tier'],'NA','NA','原808代表' if m in {a['MAG_ID'] for a in meta} else '补选阳性代表'])
for m,p,k in hits:counts[m,k].add(p)
carriers={m for m,p,k in hits};orig={x['MAG_ID'] for x in meta};taxrows=[];carrierrows=[]
def taxpart(t,p):return next((s[3:] for s in t.split(';') if s.startswith(p)), '') or 'Unclassified'
for m in carriers:
 x=md[m];tax=x['classification'];d=taxpart(tax,'d__');ph=taxpart(tax,'p__');geneids=[byko[k][1] for k in order if counts[m,k]];source='原808全局代表' if m in orig else '补选mcrA阳性代表'
 taxrows.append([m,d,ph]+[taxpart(tax,p) for p in ['c__','o__','f__','g__','s__']]+[tax,'本MAG GTDB结果' if m in orig else '同一dRep簇原代表分类代理',x['completeness'],x['contamination'],source])
 carrierrows.append([m,m+'.fasta',source,x['dRep_cluster'],int(m in sel),len(geneids),'; '.join(geneids),'; '.join(k for k in order if counts[m,k]),len({p for mm,p,k in hits if mm==m}),'直接MAG TPM可用' if m in orig else '直接MAG TPM缺失'])
taxrows.sort(key=lambda x:(x[1],x[2],x[0]));taxdict={x[0]:x for x in taxrows};carrierrows.sort(key=lambda x:(taxdict[x[0]][1],taxdict[x[0]][2],-x[5],x[0]));carrierorder=[x[0] for x in carrierrows]
save('carrying_MAGs.tsv',['MAG_ID','bin_filename','cohort','dRep_cluster','mcrA_positive_selected','target_gene_types','gene_names','KO_list','protein_records','abundance_status'],carrierrows)
save('MAG_taxonomy.tsv',['MAG_ID','domain','phylum','class','order','family','genus','species','full_GTDB_taxonomy','taxonomy_evidence','completeness_percent','contamination_percent','cohort'],taxrows)
summary=[]
for d,ph in sorted({(x[1],x[2]) for x in taxrows}):
 ms={x[0] for x in taxrows if x[1:3]==[d,ph]};summary.append([d,ph,len(ms&orig),len(ms-orig),len({md[m]['dRep_cluster'] for m in ms}),len(ms&set(sel))])
save('phylum_summary.tsv',['domain','phylum','original808_carrier_MAGs','supplement_carrier_MAGs','unique_dRep_clusters','selected_mcrA_positive_MAGs'],summary)
save('protein_KO_hits.tsv',['MAG_ID','protein_ID','KO','gene_name','evidence','score','E_value','cohort'],sorted(hits.values()))
for suffix,genes,ords in [('CH4_CO2',pathgenes,order),('reference_style',refgenes,reforder)]:
 save('all_carrier_detection_'+suffix+'.tsv',['MAG_ID','domain','phylum']+genes,[[m,taxdict[m][1],taxdict[m][2]]+[int(bool(counts[m,k])) for k in ords] for m in carrierorder])
 focusorder=sorted(sel,key=lambda m:({'ANME-2c':0,'ANME-2a':1,'ANME-1':2,'unresolved':3,'methanogen_related_unresolved':4}.get(sel[m]['clade'],5),m))
 save('mcrA_positive_detection_'+suffix+'.tsv',['MAG_ID','clade','phylum','completeness_percent']+genes,[[m,sel[m]['clade'],taxdict[m][2],md[m]['completeness']]+[int(bool(counts[m,k])) for k in ords] for m in focusorder])
save('MAG_gene_ORF_counts.tsv',['MAG_ID']+pathgenes,[[m]+[len(counts[m,k]) for k in order] for m in carrierorder])
samples=[x['sample'] for x in read(P/'tables/sample_metadata.tsv')]
community={x['KO']:x for x in read(P/'local_audit/local_community_pathway_KO_56samples.tsv')};maggenes={x['KO']:x for x in read(P/'tables/gene_TPM_56samples_mcrA_positive_MAGs.tsv')}
for name,data in [('community_gene_abundance',community),('mcrA_MAG_gene_abundance',maggenes)]:
 save(name+'.tsv',['gene_name','KO']+samples,[[byko[k][1],k]+[data.get(k,{}).get(s,'NA') for s in samples] for k in order])
save('KO_gene_mapping.tsv',['KO','gene_name','functional_group','function','CH4_CO2_order','reference_order','KEGG_source'],[[k,n,g,f,order.index(k)+1,reforder.index(k)+1,'https://www.kegg.jp/entry/'+k] for k,n,g,f in G])
shutil.copy2(P/'tables/selected_positive_MAG_TPM_56samples.tsv',T/'selected_MAG_TPM_56samples.tsv');shutil.copy2(P/'tables/sample_metadata.tsv',T/'sample_metadata.tsv')
metrics={'target_KOs':23,'original_MAGs_screened':808,'original_MAG_carriers':len(carriers&orig),'supplement_MAG_carriers':len(carriers-orig),'carrier_records_total':len(carriers),'unique_carrier_dRep_clusters':len({md[m]['dRep_cluster'] for m in carriers}),'mcrA_positive_representatives':16,'taxonomic_phyla':len({x[2] for x in taxrows}),'protein_KO_records':len(hits)}
(O/'summary.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2));(O/'plot_config.json').write_text(json.dumps({'genes':pathgenes,'KO_order':order,'reference_genes':refgenes,'reference_KO_order':reforder,'focus_order':focusorder},ensure_ascii=False))
print(json.dumps(metrics,ensure_ascii=False))
