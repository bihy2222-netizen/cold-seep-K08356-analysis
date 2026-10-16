"""Verify the archived cold-seep methane metabolism summary; Python standard library only."""
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def read(name):
    with (ROOT / 'tables' / name).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

carriers = read('carrying_MAGs.tsv')
taxonomy = read('MAG_taxonomy.tsv')
focus = read('mcrA_positive_detection_CH4_CO2.tsv')
mapping = read('KO_gene_mapping.tsv')
hits = read('protein_KO_hits.tsv')
samples = read('sample_metadata.tsv')
original = {r['MAG_ID'] for r in read('cold_seep_all808_taxonomy_quality.tsv')}
ids = {r['MAG_ID'] for r in carriers}
assert len(original) == 808
assert len(carriers) == len(ids) == 347
assert len(ids & original) == 344 and len(ids - original) == 3
assert len({r['dRep_cluster'] for r in carriers}) == 344
assert ids == {r['MAG_ID'] for r in taxonomy}
assert len({r['phylum'] for r in taxonomy}) == 43
assert len(focus) == 16 and all(r['phylum'] == 'Halobacteriota' for r in focus)
assert Counter(r['clade'] for r in focus) == Counter({'ANME-1': 2, 'ANME-2a': 2, 'ANME-2c': 2, 'unresolved': 8, 'methanogen_related_unresolved': 2})
assert len(mapping) == len({r['KO'] for r in mapping}) == 23
assert 'K00204' not in {r['KO'] for r in mapping}
assert len(hits) == len({(r['MAG_ID'], r['protein_ID'], r['KO']) for r in hits}) == 1876
assert Counter(r['habitat'] for r in samples) == Counter({'NS': 8, 'IS': 21, 'AS': 14, 'ES': 13})
for name in ['community_gene_abundance.tsv', 'mcrA_MAG_gene_abundance.tsv']:
    rows = read(name)
    assert len(rows) == 23
    assert list(rows[0])[2:] == [r['sample'] for r in samples]
for name, metadata in [('all_carrier_detection_CH4_CO2.tsv', 3), ('all_carrier_detection_reference_style.tsv', 3), ('mcrA_positive_detection_CH4_CO2.tsv', 4), ('mcrA_positive_detection_reference_style.tsv', 4)]:
    for row in read(name):
        values = list(row.values())[metadata:]
        assert len(values) == 23 and set(values) <= {'0', '1'} and '1' in values
result = {'status': 'PASS', 'MAG_records': 347, 'independent_dRep_clusters': 344, 'mcrA_positive_representatives': 16, 'target_KOs': 23, 'samples': 56, 'protein_KO_records': 1876}
(ROOT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result))
