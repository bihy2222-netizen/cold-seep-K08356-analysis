#!/usr/bin/env python3
"""Separate true junction-spanning reads from paired-end boundary support."""

import argparse,csv,re,subprocess
from collections import defaultdict

CIGAR=re.compile(r"(\d+)([MIDNSHP=X])")
def reference_end(start,cigar):return start+sum(int(n) for n,op in CIGAR.findall(cigar) if op in "MDN=X")-1
def original_contig(ref):
    if "|original=" in ref:return ref.split("|original=",1)[1].split("|",1)[0]
    if "|contig=" in ref:return ref.split("|contig=",1)[1].split("|",1)[0]
    return ref
def overlaps(read,start,end):return read["start"]<=end and read["end"]>=start

def main():
    p=argparse.ArgumentParser();p.add_argument("--bam",required=True);p.add_argument("--boundaries",required=True)
    p.add_argument("--crosswalk");p.add_argument("--out",required=True);a=p.parse_args()
    boundaries=list(csv.DictReader(open(a.boundaries),delimiter="\t"));wanted={r["contig_id"] for r in boundaries}
    reference_to_original={}
    if a.crosswalk:
        with open(a.crosswalk,newline="") as handle:
            reference_to_original={r["reference_id"]:r["original_contig_id"] for r in csv.DictReader(handle,delimiter="\t")}
    reads=defaultdict(dict);excluded={"secondary":0,"supplementary":0}
    proc=subprocess.Popen(["samtools","view",a.bam],stdout=subprocess.PIPE,text=True)
    for line in proc.stdout:
        f=line.rstrip().split("\t");flag=int(f[1]);contig=reference_to_original.get(f[2],original_contig(f[2]))
        if flag&256:excluded["secondary"]+=1;continue
        if flag&2048:excluded["supplementary"]+=1;continue
        if flag&4 or contig not in wanted or f[5]=="*":continue
        mate=1 if flag&64 else 2 if flag&128 else 0
        reads[(contig,f[0])][mate]={"start":int(f[3]),"end":reference_end(int(f[3]),f[5]),"mapq":int(f[4]),
            "flag":flag,"tlen":abs(int(f[8])),"reverse":bool(flag&16),"duplicate":bool(flag&1024)}
    if proc.wait()!=0:raise SystemExit("samtools view failed")
    output=[]
    for b in boundaries:
        contig=b["contig_id"];gap_left=int(b["genomic_left_gene_end"]);gap_right=int(b["genomic_right_gene_start"])
        first=(int(b["first_gene_start"]),int(b["first_gene_end"]));second=(int(b["second_gene_start"]),int(b["second_gene_end"]))
        c=defaultdict(int);inserts=[]
        for (read_contig,_),mates in reads.items():
            if read_contig!=contig:continue
            spanning=[m for m in mates.values() if m["start"]<=gap_left and m["end"]>=gap_right]
            if spanning:
                c["spanning_read_count"]+=1
                if any(m["mapq"]>=20 for m in spanning):c["MAPQ20_spanning_reads"]+=1
                if all(not m["duplicate"] for m in spanning):c["duplicate_excluded_spanning_reads"]+=1
            if 1 in mates and 2 in mates:
                one,two=mates[1],mates[2]
                pair_spans=(overlaps(one,*first) and overlaps(two,*second)) or (overlaps(two,*first) and overlaps(one,*second))
                if pair_spans:
                    c["spanning_pair_count"]+=1;inserts.append(max(one["tlen"],two["tlen"]))
                    if min(one["mapq"],two["mapq"])>=20:c["MAPQ20_spanning_pairs"]+=1
                    if one["flag"]&2 and two["flag"]&2:c["proper_pair_count"]+=1
                    left_read,right_read=sorted((one,two),key=lambda read:read["start"])
                    if not left_read["reverse"] and right_read["reverse"]:c["orientation_consistent_count"]+=1
                    if not one["duplicate"] and not two["duplicate"]:c["duplicate_excluded_spanning_pairs"]+=1
        output.append({"boundary_id":f"{b['candidate_id']}|{b['first_role_transcriptional']}--{b['second_role_transcriptional']}",
            "candidate_id":b["candidate_id"],"first_role":b["first_role_transcriptional"],"second_role":b["second_role_transcriptional"],
            "contig_id":contig,"strand":b["strand"],"intergenic_distance":b["intergenic_bp"],
            **{k:c[k] for k in ("spanning_read_count","spanning_pair_count","MAPQ20_spanning_reads","MAPQ20_spanning_pairs",
                "proper_pair_count","orientation_consistent_count","duplicate_excluded_spanning_reads","duplicate_excluded_spanning_pairs")},
            "median_insert_size":sorted(inserts)[len(inserts)//2] if inserts else "",
            "secondary_excluded_count_global":excluded["secondary"],"supplementary_excluded_count_global":excluded["supplementary"],
            "evidence_label":"DNA_presence_control_not_cotranscription"})
    fields=list(output[0]) if output else []
    with open(a.out,"w",newline="") as h:w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(output)
if __name__=="__main__":main()
