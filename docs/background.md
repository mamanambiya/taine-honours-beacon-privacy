# Background reading

The proposal's references, grouped by the job each one does in the argument.
Read the four marked **core** first — the rest support the motivation and can be
read as the write-up needs them.

## The attack this project reimplements

- **core** — Shringarpure, S. S., & Bustamante, C. D. (2015). Privacy Risks from
  Genomic Data-Sharing Beacons. *Am J Hum Genet*, 97(5), 631–646.
  [doi:10.1016/j.ajhg.2015.09.010](https://doi.org/10.1016/j.ajhg.2015.09.010)
  — the LRT membership-inference attack. The whole method section derives from
  this. Get the exact statistic and the required-`n` expression from here.

- Ayoz, K., Ayday, E., & Cicek, A. E. (2021). Genome Reconstruction Attacks
  Against Genomic Data-Sharing Beacons. *PoPETs*, 2021(3), 28–48.
  [doi:10.2478/popets-2021-0036](https://doi.org/10.2478/popets-2021-0036)
  — exploits LD correlations; the reason a defence tuned only against
  membership inference is insufficient.

- Saleem, K., Cicek, A. E., & Sav, S. (2025). Beacon Reconstruction Attack.
  *Bioinformatics*, 41(6).
  [doi:10.1093/bioinformatics/btaf273](https://doi.org/10.1093/bioinformatics/btaf273)
  — the current state of the art; optimisation-based reconstruction from
  summary statistics.

## Proposed defences, and why the proposal finds them unsatisfying

- Aziz, M. M. Al, et al. (2017). Aftermath of Bustamante attack on genomic
  beacon service. *BMC Med Genomics*, 10(S2), 43.
  [doi:10.1186/s12920-017-0278-x](https://doi.org/10.1186/s12920-017-0278-x)
  — risk-analysis threshold plus deliberately inaccurate results.

- Raisaro, J. L., et al. (2017). Addressing Beacon re-identification attacks.
  *JAMIA*, 24(4), 799–805.
  [doi:10.1093/jamia/ocw167](https://doi.org/10.1093/jamia/ocw167)
  — the query-budget defence.

- Dyke, S. O., Dove, E. S., & Knoppers, B. M. (2016). Sharing health-related
  data: a privacy test? *npj Genomic Medicine*, 1(1), 16024.
  [doi:10.1038/npjgenmed.2016.24](https://doi.org/10.1038/npjgenmed.2016.24)
  — GA4GH's contextual Data Sharing Privacy Test.

Both technical defences trade utility for privacy — inaccurate answers, or no
answers past a budget. The proposal's argument is that the size of that trade
has never been costed for African cohorts, because the attack power was never
measured on one.

## Beacon as a protocol

- **core** — Rambla, J., et al. (2022). Beacon v2 and Beacon networks: a
  "lingua franca" for federated data discovery. *Human Mutation*.
  [doi:10.1002/humu.24369](https://doi.org/10.1002/humu.24369)
  — the v2 spec this deployment implements; granularity levels and access tiers.

- Fiume, M., et al. (2019). Federated discovery and sharing of genomic data
  using Beacons. *Nat Biotechnol*, 37(3), 220–224.
  [doi:10.1038/s41587-019-0046-x](https://doi.org/10.1038/s41587-019-0046-x)

- Rueda, M., Ariosa, R., Moldes, M., & Rambla, J. (2022). Beacon v2 Reference
  Implementation. *Bioinformatics*, 38(19), 4656–4657.
  [doi:10.1093/bioinformatics/btac568](https://doi.org/10.1093/bioinformatics/btac568)

- Naro, D., Delgado, J., Llorente, S., & Hernaez, M. (2024). Use of Beacon v2
  for Improving Genomics Based Research in a Clinical Setting.
  [doi:10.3233/SHTI240636](https://doi.org/10.3233/SHTI240636)

- Knoppers, B. M. (2014). Framework for responsible sharing of genomic and
  health-related data. *The HUGO Journal*, 8(1), 3.
  [doi:10.1186/s11568-014-0003-1](https://doi.org/10.1186/s11568-014-0003-1)

## Why the African context changes the question

- **core** — Tiffin, N. (2019). Potential risks and solutions for sharing genome
  summary data from African populations. *BMC Med Genomics*, 12(1), 152.
  [doi:10.1186/s12920-019-0604-6](https://doi.org/10.1186/s12920-019-0604-6)
  — the closest existing work to this project's question. Read early; it frames
  the harm side of the argument.

- **core** — Sirugo, G., Williams, S. M., & Tishkoff, S. A. (2019). The Missing
  Diversity in Human Genetic Studies. *Cell*, 177(1), 26–31.
  [doi:10.1016/j.cell.2019.02.048](https://doi.org/10.1016/j.cell.2019.02.048)

- Ojewunmi, O. O., & Fatumo, S. (2025). Driving global health equity and
  precision medicine through African genomic data. *Hum Mol Genet*.
  [doi:10.1093/hmg/ddaf025](https://doi.org/10.1093/hmg/ddaf025)

- Omotoso, O. E., et al. (2022). Bridging the genomic data gap in Africa.
  *Globalization and Health*, 18(1), 103.
  [doi:10.1186/s12992-022-00898-2](https://doi.org/10.1186/s12992-022-00898-2)

- Lumaka, A., et al. (2022). Increasing African genomic data generation and
  sharing. *Orphanet J Rare Dis*, 17(1), 230.
  [doi:10.1186/s13023-022-02391-w](https://doi.org/10.1186/s13023-022-02391-w)

- Fatumo, S., & Inouye, M. (2023). African genomes hold the key to accurate
  genetic risk prediction. *Nat Hum Behav*, 7(3), 295–296.
  [doi:10.1038/s41562-023-01549-1](https://doi.org/10.1038/s41562-023-01549-1)

- Kamp, M., Krause, A., & Ramsay, M. (2021). Has translational genomics come of
  age in Africa? *Hum Mol Genet*, 30(R2), R164–R173.
  [doi:10.1093/hmg/ddab180](https://doi.org/10.1093/hmg/ddab180)

- Madden, E. B., et al. (2024). Advancing genomics to improve health equity.
  *Nat Genet*, 56(5), 752–757.
  [doi:10.1038/s41588-024-01711-z](https://doi.org/10.1038/s41588-024-01711-z)

- Mills, M. C., & Rahal, C. (2019). A scientometric review of GWAS.
  *Commun Biol*, 2(1), 9.
  [doi:10.1038/s42003-018-0261-x](https://doi.org/10.1038/s42003-018-0261-x)

- Manrai, A. K., et al. (2016). Genetic Misdiagnoses and the Potential for
  Health Disparities. *NEJM*, 375(7), 655–665.
  [doi:10.1056/NEJMsa1507092](https://doi.org/10.1056/NEJMsa1507092)

## Gaps worth checking before the literature review is written

The proposal's citation set stops at the attack/defence literature and the
African-genomics motivation. Two areas are not represented and a reviewer will
look for them:

1. **Beacon-specific differential privacy** — work post-2021 applying DP or
   noise calibration to beacon responses, as an alternative to the query-budget
   and censoring defences already cited.
2. **Empirical membership inference on non-European cohorts** — the proposal
   asserts this has not been done for African data. That is the study's
   novelty claim, so it needs a documented search rather than an assumption.
   Establish it with a recorded search strategy; "I did not find one" and
   "there is not one" are different sentences, and only the second supports a
   novelty claim.
