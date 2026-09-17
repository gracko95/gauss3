# PRS-CS reference-panel mixing scripts

`PRScs_mix.py` and `v2parse_genet_mix.py` extend [PRS-CS](https://github.com/getian107/PRScs)
(Ge et al. 2019) to compute posterior SNP effects from a weighted blend of multiple
reference LD panels (e.g. EUR + ASN), instead of a single population panel. They are
not part of upstream PRS-CS and were not previously under version control anywhere.

To run, place these alongside an unmodified PRS-CS checkout so `mcmc_gtb.py` and
`gigrnd.py` are importable (`PRScs_mix.py` imports `v2parse_genet_mix`, `mcmc_gtb`,
`gigrnd`).

`PRScs_mix.py`'s `main()` intersects every `--ref_dir` panel's SNP universe up front
before building the SNP list. Earlier versions built the SNP list from only the
first-listed panel, so a SNP missing from a smaller panel (e.g. ASN, which covers
fewer chr22 SNPs than EUR due to population-specific MAF filtering) was silently
dropped inside `parse_ldblk_mix`'s per-block LD computation instead of being excluded
from the SNP list -- leaving it with a posterior beta of exactly 0.0. Because the
sampler's shared noise-scale parameter in `mcmc_gtb.py` is drawn using the total SNP
count, those phantom zeroed SNPs also biased every other SNP's estimate.
