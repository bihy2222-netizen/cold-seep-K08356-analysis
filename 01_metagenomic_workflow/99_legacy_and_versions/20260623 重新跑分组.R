##############20260623师兄重新跑了derep#################
重命名bin名称
(base) ps@ps:~/ps1/data/bihongyu/cold_seep/illu/binning$ nohup bash rename_bins.sh &
  去冗余
(base) ps@ps:~/ps1/data/bihongyu/cold_seep/illu/binning/all_50_10_dRep$ nohup dRep dereplicate ./all_50_10_dRep -g ./all_refined_bin/*.fasta -p 150 --ignoreGenomeQuality &