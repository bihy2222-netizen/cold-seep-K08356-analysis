# Standard METABOLIC-G workflow

The reproducible K08356 MAG workflow is documented in:

`01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/README.md`

The key distinction is analytical scope:

- `all_50_10_dRep/.../bin_METABOLIC` is a completed 808-MAG background run.
- The original subset has 48 K08356 proteins from 47 unique MAGs.
- Habitat-wise dereplication yields 49 K08356 proteins from 48 unique MAGs and
  adds `R2111_N500_0-10_bin13`, so a separate run is required.

Use the K08356 hit table to select unique MAG-level protein inputs. Do not pass
the 49 individual K08356 proteins to METABOLIC-G; the tool requires each MAG's
complete protein set (`*.cds.faa`) to reconstruct metabolic modules and element
cycles. `S1_9-12_bin1` has two K08356 proteins and is included once, giving 48
MAG inputs.

The 48-MAG METABOLIC-G v4.0 run completed on 2026-08-01 in 35 minutes 28
seconds. All cardinality and completion checks passed; see the workflow
manifest directory for the exact inputs, parameters, and validation summary.
