#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict

def read(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def family(s):
    if s.startswith(("IdrP_like","IdrP1","IdrP2")): return "P_like"
    if s.startswith("IdrB"): return "IdrB_related"
    if s.startswith("AioB"): return "canonical_AioB"
    return "other"

def accepted(fam, e, bits, qcov, scov):
    if fam in {"IdrB_related","canonical_AioB"}:
        return e <= 1e-5 and bits >= 50 and qcov >= 50 and scov >= 35
    if fam == "P_like":
        return e <= 1e-10 and bits >= 80 and qcov >= 45 and scov >= 35
    return False

ap=argparse.ArgumentParser()
ap.add_argument("--neighborhood-long",required=True)
ap.add_argument("--neighborhood-summary",required=True)
ap.add_argument("--diamond",required=True)
ap.add_argument("--out-hits",required=True)
ap.add_argument("--out-summary",required=True)
a=ap.parse_args()

long=read(a.neighborhood_long); status={r["target_id"]:r for r in read(a.neighborhood_summary)}
meta={r["neighbor_raw_id"]:r for r in long}
hits=[]
with open(a.diamond,newline="") as h:
    for p in csv.reader(h,delimiter="\t"):
        if len(p)<14 or p[0] not in meta: continue
        q,s=p[0],p[1]; alen=float(p[3]); qlen=float(p[12]); slen=float(p[13])
        fam=family(s); e=float(p[10]); bits=float(p[11])
        qcov=100*alen/qlen if qlen else 0; scov=100*alen/slen if slen else 0
        r=dict(meta[q],qseqid=q,sseqid=s,ref_family=fam,evalue=e,bitscore=bits,
               qcov_pct=round(qcov,2),scov_pct=round(scov,2))
        r["accepted"]="yes" if accepted(fam,e,bits,qcov,scov) else "no"
        hits.append(r)

hf=list(hits[0].keys()) if hits else ["target_id","qseqid","sseqid","ref_family","accepted"]
with open(a.out_hits,"w",newline="") as h:
    w=csv.DictWriter(h,delimiter="\t",fieldnames=hf,extrasaction="ignore"); w.writeheader(); w.writerows(hits)

by=defaultdict(list)
for r in hits:
    if r["accepted"]=="yes": by[r["target_id"]].append(r)

out=[]
for tid, st in sorted(status.items()):
    acc=by.get(tid,[])
    b=sorted([r for r in acc if r["ref_family"] in {"IdrB_related","canonical_AioB"}],
             key=lambda r:-r["bitscore"])
    pbest={}
    for r in [x for x in acc if x["ref_family"]=="P_like"]:
        if r["qseqid"] not in pbest or r["bitscore"]>pbest[r["qseqid"]]["bitscore"]:
            pbest[r["qseqid"]]=r
    ps=sorted(pbest.values(),key=lambda r:int(r["relative_position"]))
    support=([b[0]] if b else [])+ps[:2]
    positions=[0]+[int(r["relative_position"]) for r in support]
    compact=(len(support)>=3 and max(positions)-min(positions)<=3)
    strands=[r["strand"] for r in support]
    cotranscribed=(len(support)>=3 and len(set(strands))==1 and compact)
    observable=st.get("edge_censoring")=="complete_observable"
    points=[bool(b),len(ps)>=1,len(ps)>=2,compact,observable,cotranscribed]
    raw=sum(points); available=6
    if not observable:
        syn=st.get("edge_censoring","edge_censored")
    elif b and len(ps)>=2 and compact:
        syn="complete"
    elif b or ps:
        syn="partial_but_observable"
    else:
        syn="absent"
    out.append(dict(target_id=tid,group=st.get("group",""),genome_id=st.get("genome_id",""),
        edge_censoring=st.get("edge_censoring",""),B_related_present="yes" if b else "no",
        B_reference=b[0]["sseqid"] if b else "",P_like_count=len(ps),
        P_like_references=";".join(r["sseqid"] for r in ps),
        compact_A_B_P_P_order="yes" if compact else "no",
        direction_and_spacing_support_cotranscription="yes" if cotranscribed else "no",
        signal_localization_score="NA_not_run",synteny_score_raw=raw,
        synteny_score_available=available,synteny_score_normalized=round(raw/available,3),
        synteny_status=syn))
fields=list(out[0].keys()) if out else ["target_id","synteny_status"]
with open(a.out_summary,"w",newline="") as h:
    w=csv.DictWriter(h,delimiter="\t",fieldnames=fields);w.writeheader();w.writerows(out)
