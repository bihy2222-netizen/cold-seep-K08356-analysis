from pathlib import Path
import json,math
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
O=Path(__file__).resolve().parent;T=O/'tables';F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'hatch.linewidth':.45})
cfg=json.loads((O/'plot_config.json').read_text());clcolors={'ANME-2c':'#A8CE9A','ANME-2a':'#8FA8D0','ANME-1':'#F4AC82','unresolved':'#8798AA','methanogen_related_unresolved':'#8798AA'}
def save(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(F/(name+'.'+ext),dpi=200,bbox_inches='tight')
 plt.close(fig)
def pathway(ax):
 ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
 nodes=['CH₄','CH₃–S–CoM','methyl-H₄MPT','methylene-H₄MPT','methenyl-H₄MPT','formyl-H₄MPT','formyl-MFR','CO₂']
 xx=np.linspace(.025,.975,8)
 for x,txt in zip(xx,nodes):ax.text(x,.65,txt,ha='center',va='center',fontsize=9)
 for i,label in enumerate(['Mcr','Mtr','Mer','Mtd','Mch','Ftr','Fwd/Fmd']):
  a,b=xx[i]+.042,xx[i+1]-.042;ax.annotate('',xy=(b,.64),xytext=(a,.64),arrowprops={'arrowstyle':'->','color':'#354A5E','lw':1.2});ax.text((xx[i]+xx[i+1])/2,.35,label,ha='center',fontsize=9,color='#354A5E')
 ax.text(.5,.03,'Conceptual CH₄ → CO₂ carbon pathway; HdrA2/B2/C2 shown separately as an electron-transfer/cofactor-recycling group.',ha='center',fontsize=8,color='#555555')
def focus(style):
 data=pd.read_csv(T/f'mcrA_positive_detection_{style}.tsv',sep='\t');genes=cfg['reference_genes'] if style=='reference_style' else cfg['genes'];n=len(data)
 fig=plt.figure(figsize=(16.8,8.7));gs=fig.add_gridspec(2,1,height_ratios=[1,5.4],hspace=.22);pathway(fig.add_subplot(gs[0]));ax=fig.add_subplot(gs[1])
 ax.set_xlim(-11.4,23);ax.set_ylim(n+.2,-3.4);ax.axis('off')
 groups=[('Mcr',3),('Mtr',8),('Mer',1),('Mtd',1),('Mch',1),('Ftr',1),('Fwd/Fmd',5),('Hdr2',3)];offset=0
 for label,width in groups:
  ax.add_patch(Rectangle((offset,-3.2),width,.58,facecolor='#EDF1F5',edgecolor='#BAC4CC',lw=.5));ax.text(offset+width/2,-2.91,label,ha='center',va='center',fontsize=9,fontweight='bold');offset+=width
 for j,g in enumerate(genes):ax.text(j+.5,-.28,g,rotation=90,ha='left',va='bottom',fontsize=8)
 ax.text(-11.3,-.35,'mcrA-positive representative MAG',va='bottom',fontweight='bold',fontsize=9);ax.text(-1.03,-.35,'Comp.\n(%)',ha='center',va='bottom',fontsize=8)
 prev=None
 for i,(_,r) in enumerate(data.iterrows()):
  color=clcolors.get(r['clade'],'#8798AA')
  if prev is not None and prev!=r['clade']:ax.axhline(i-.06,color='white',lw=2)
  prev=r['clade'];ax.add_patch(Rectangle((-11.35,i),9.65,1,facecolor=color,alpha=.85,edgecolor='white',lw=.6));ax.text(-11.16,i+.5,r.MAG_ID,ha='left',va='center',fontsize=8.4)
  c=r['completeness_percent'];ct='NA' if pd.isna(c) else f'{float(c):.0f}';ax.add_patch(Rectangle((-1.48,i),.87,1,facecolor=color,edgecolor='white',lw=.6));ax.text(-1.045,i+.5,ct,ha='center',va='center',fontsize=8)
  for j,g in enumerate(genes):
   detected=int(r[g])==1;ax.add_patch(Rectangle((j,i),1,1,facecolor=color if detected else '#E1E1E1',edgecolor='#9A9A9A',lw=.45,hatch=None if detected else '///'))
 handles=[Patch(facecolor=clcolors[g],label=g) for g in ['ANME-2c','ANME-2a','ANME-1']]+[Patch(facecolor='#8798AA',label='Unresolved McrA group'),Patch(facecolor='#E1E1E1',edgecolor='#999999',hatch='///',label='Not detected')]
 fig.legend(handles=handles,ncol=5,loc='lower center',bbox_to_anchor=(.52,.055),frameon=False,fontsize=8)
 fig.suptitle(('A. Reference-style MAG gene-detection matrix (no tree)' if style=='reference_style' else 'B. MAG gene detection in CH₄ → CO₂ pathway order (no tree)'),x=.06,ha='left',fontsize=13)
 fig.text(.06,.023,'Filled cells: genomic detection, not abundance. All 16 hosts: GTDB Halobacteriota; 3 supplemental MAGs use representative-proxy taxonomy. No operon boundaries inferred.',fontsize=8)
 save(fig,'01_reference_style_MAG_gene_names' if style=='reference_style' else '02_CH4_CO2_MAG_gene_names')
def abundance(filename,title,name,source):
 a=pd.read_csv(T/filename,sep='\t');meta=pd.read_csv(T/'sample_metadata.tsv',sep='\t');ss=meta['sample'].tolist();values=a[ss].to_numpy(float)
 fig,ax=plt.subplots(figsize=(17,9));fig.subplots_adjust(bottom=.22,top=.89,left=.14,right=.95);cmap=plt.get_cmap('magma').copy();cmap.set_bad('#CBD5E1');im=ax.imshow(np.ma.masked_invalid(np.log10(1+values)),aspect='auto',cmap=cmap,vmin=0)
 ax.set_yticks(range(len(a)),a.gene_name,fontsize=9);ax.set_xticks(range(56),ss,rotation=90,fontsize=7);ax.set_title(title,loc='left',pad=32,fontsize=12)
 for g,c in {'NS':'#3B7EA1','IS':'#E7AD3E','AS':'#C65252','ES':'#578A6A'}.items():
  ix=np.flatnonzero(meta.habitat.to_numpy()==g);ax.plot([min(ix)-.5,max(ix)+.5],[-1,-1],lw=5,color=c,clip_on=False);ax.text(ix.mean(),-1.65,g,ha='center',clip_on=False)
  if min(ix)>0:ax.axvline(min(ix)-.5,color='white',lw=.8)
 # Reaction groups follow CH4 oxidation order, rather than KO number.
 breaks=[(0,3,'Mcr'),(3,11,'Mtr'),(11,12,'Mer'),(12,13,'Mtd'),(13,14,'Mch'),(14,15,'Ftr'),(15,20,'Fwd/Fmd'),(20,23,'Hdr2')]
 for start,end,label in breaks:
  if start:ax.axhline(start-.5,color='white',lw=.9)
 fig.colorbar(im,ax=ax,pad=.015,shrink=.7,label='log10('+('source-reported KO abundance' if source=='community' else 'gene TPM')+' + 1)')
 fig.text(.14,.025,('Community KO aggregation; 56 samples; not attributed to individual MAGs.' if source=='community' else 'Gene TPM contributions from 13 directly quantified mcrA-positive representatives; 3 supplemental MAGs not included.')+' Grey = unavailable quantification.',fontsize=8)
 save(fig,name)
def all_carriers():
 a=pd.read_csv(T/'all_carrier_detection_CH4_CO2.tsv',sep='\t');genes=cfg['genes'];out=F/'05_all_347_carrier_MAGs_no_tree_paginated.pdf'
 with PdfPages(out) as pdf:
  for begin in range(0,len(a),40):
   sub=a.iloc[begin:begin+40];fig,ax=plt.subplots(figsize=(17,12));fig.subplots_adjust(left=.37,bottom=.18,top=.94,right=.99)
   vals=sub[genes].to_numpy();ax.imshow(vals,aspect='auto',vmin=0,vmax=1,cmap=ListedColormap(['#E4E4E4','#527C9E']))
   ax.set_yticks(range(len(sub)),[f'{m} | {p}' for m,p in zip(sub.MAG_ID,sub.phylum)],fontsize=8);ax.set_xticks(range(23),genes,rotation=90,fontsize=8);ax.set_title(f'All target-gene carriers: {begin+1}–{begin+len(sub)} of {len(a)} (gene detection, no tree)',fontsize=11,loc='left')
   for x in np.arange(-.5,23,1):ax.axvline(x,color='white',lw=.4)
   fig.text(.37,.035,'Blue = detected; grey = not detected. Any of 23 genes qualifies as a carrier; this does not identify an ANME or a complete methane pathway.',fontsize=8)
   pdf.savefig(fig);plt.close(fig)
focus('reference_style');focus('CH4_CO2')
abundance('community_gene_abundance.tsv','C. Community gene-associated KO abundance in CH₄ → CO₂ pathway order','03_CH4_CO2_community_gene_names','community')
abundance('mcrA_MAG_gene_abundance.tsv','D. mcrA-positive MAG gene TPM in CH₄ → CO₂ pathway order','04_CH4_CO2_mcrA_MAG_gene_names','MAG')
all_carriers();print('FIGURES_COMPLETE')
