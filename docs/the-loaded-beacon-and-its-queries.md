# The loaded beacon, and every query that explores it

What the beacon on `afrigend-all-dev` is, what it holds, and the queries that
work against it. Every command below was run on 2026-10-28 unless it says
"from code, not run". The beacon is under active development; re-run before
relying on a number.

How to reach it, check it and attack it is in
[`querying-and-attacking-the-loaded-beacons.md`](querying-and-attacking-the-loaded-beacons.md).
Its operational limits (rate limit, GRCh38 only, N is public) are in
[`beacon-constraints.md`](beacon-constraints.md).

All commands run on this Mac through the ssh alias `afrigend-all-dev`; the API
itself listens on `127.0.0.1:8086` on the VM only. Each prints JSON through
`python3 -m json.tool`.

## 1. What it is

| | |
| --- | --- |
| Beacon id | `org.afrigen-d.beacon`, GA4GH Beacon v2.0.0 |
| Mode | boolean by default (`"version": "2.0.0-boolean"`), public, no login |
| Assembly | GRCh38 only; GRCh37 answers 501 |
| Datasets | 12: `<POP>_M65_S1` for CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI and the pooled AFR, EUR |
| Variants | 9,859,803 across the 12 (chromosomes 20 to 22; each record means a member carries the alternate allele) |
| Per dataset | 65 members of the 1000 Genomes NYGC 30x panel; 20 held-out controls are not loaded |
| Individuals, biosamples, analyses | none: `/api/individuals` answers `exists: false`, 0 results |
| Code | `AfriGen-D/variant-checker-beacon` (`beacon_api/`, boolean views in `views_boolean.py`) |

## 2. Endpoints

From `beacon_api/urls_boolean.py`. "Run" means the probe below was run.

| Path | What it gives | Run |
| --- | --- | --- |
| `/api/` , `/api/info`, `/api/service-info` | beacon identity and the dataset list | yes |
| `/api/map` | endpoint map | yes |
| `/api/configuration` | entry types and security level (`PUBLIC`, default granularity `boolean`) | yes |
| `/api/entry_types` | the entry types served | yes |
| `/api/health` | database and cache status | yes |
| `/api/datasets` | all 12 datasets with sizes | yes |
| `/api/datasets/<id>` | one dataset | yes |
| `/api/datasets/<id>/<entry_type>` | a query scoped to one dataset | from code, not run |
| `/api/filtering_terms` | filtering terms (none defined) | yes |
| `/api/g_variants`, `/api/query/variants` | the variant query, the main one | yes |
| `/api/individuals`, `/api/query/individuals` | individual query (no individuals loaded) | yes |
| `/api/biosamples`, `/api/analyses`, `/api/runs`, `/api/cohorts` | spec-shaped empty lists | from code, not run |

## 3. Describe the beacon

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/info'" | python3 -m json.tool
```

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/configuration'" | python3 -m json.tool
```

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/entry_types'" | python3 -m json.tool
```

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/map'" | python3 -m json.tool
```

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/health'" | python3 -m json.tool
```

Expect `"status": "healthy"` from the last one.

## 4. List the datasets

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/datasets'" | python3 -m json.tool
```

Expect `numTotalResults: 12`, one entry per `<POP>_M65_S1`.

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/datasets/CEU_M65_S1'" | python3 -m json.tool
```

Expect one result with the dataset's name and description.

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/filtering_terms'" | python3 -m json.tool
```

Expect `numTotalResults: 0`: no filtering terms are defined, so there is
nothing to filter on beyond the variant parameters below.

## 5. Query variants

Parameters of `/api/g_variants`, from `views_boolean.py` and
`query_vocabulary.py`:

| Parameter | Meaning |
| --- | --- |
| `assemblyId` | `GRCh38`; anything else answers 501 |
| `referenceName` | chromosome without `chr`; required whenever a position is given (400 otherwise, from code) |
| `start` | 0-based, so VCF position minus 1 |
| `end` | with `start`, a range; a range needs both |
| `referenceBases`, `alternateBases` | the alleles |
| `variantType` | `SNV` (`SNP` is accepted and means `SNV`), also `DEL INS DUP INV CNV`; the loaded data holds only SNVs |
| `datasetIds` | restrict the query to one dataset (see section 6) |
| `requestedGranularity` | `boolean` (default), `aggregated`, `record` |

**One allele in one dataset.** The standard question:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709680&referenceBases=G&alternateBases=C&datasetIds=CEU_M65_S1'" | python3 -m json.tool
```

Expect `responseSummary.exists: true`, `numTotalResults: 1`.

**A range.** Does any variant sit between two positions, in one dataset:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709600&end=709800&datasetIds=CEU_M65_S1'" | python3 -m json.tool
```

Expect `exists: true` (it includes position 709681).

**A range by type.** The same range, SNVs only:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709600&end=709800&variantType=SNV&datasetIds=CEU_M65_S1'" | python3 -m json.tool
```

**A position without alleles, over all datasets.** No `datasetIds`, no alleles:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709680'" | python3 -m json.tool
```

Expect `numTotalResults: 12`: all twelve datasets hold a variant starting there.

**Which datasets hold an allele.** Drop `datasetIds` and read the `resultSets` that
are `exists: true`:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709680&referenceBases=G&alternateBases=C'" | python3 -c "import sys,json;d=json.load(sys.stdin);print(d['responseSummary']);print([r['id'] for r in d['response']['resultSets'] if r['exists']])"
```

Expect `numTotalResults: 12` and all twelve ids for this variant.

**Everything.** No parameters is the spec's "all entries" query:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants'" | python3 -m json.tool
```

Expect `numTotalResults: 9859803`, the total over the 12 datasets, and an empty
`resultSets`. It is answered from precomputed sizes, not by scanning.

**Allele frequency.** `requestedGranularity=aggregated` adds a frequency:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709680&referenceBases=G&alternateBases=C&datasetIds=CEU_M65_S1&requestedGranularity=aggregated'" | python3 -m json.tool
```

Expect `returnedGranularity: "aggregated"` and, in `resultSets[0].results`,
`frequencyInPopulations[].frequencies[].alleleFrequency` (it was `0.1` for this
variant in CEU). The code rounds the frequency and suppresses small cells
before publishing (`beacon_api/privacy.py`), because an exact frequency
reveals the carrier count. Note for the study: this dev beacon answers more than
yes or no when asked for `aggregated`; the live production beacon was
measured boolean-only on 2026-08-21 (`beacon-constraints.md`).

**A build it does not hold:**

```bash
ssh afrigend-all-dev "curl -s -w '\nHTTP %{http_code}\n' 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh37&referenceName=20&start=709680&referenceBases=G&alternateBases=C&datasetIds=CEU_M65_S1'"
```

Expect `HTTP 501` and "This beacon holds no data for assembly GRCh37".

**Individuals.** None are loaded:

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/individuals'" | python3 -m json.tool
```

Expect `exists: false`, `numTotalResults: 0`.

## 6. Reading the response

`datasetIds` restricts the query to the dataset you name, but the response
**always lists all twelve datasets** under `response.resultSets`. Only the named
dataset is evaluated; the others show `exists: false` and `resultsCount: 0`
because they were not asked, not because they lack the variant. Measured with one
variant:

| Query | `numTotalResults` | Datasets with `exists: true` |
| --- | --- | --- |
| no `datasetIds` | 12 | all 12 |
| `datasetIds=CEU_M65_S1` | 1 | CEU only |
| `datasetIds=ESN_M65_S1` | 1 | ESN only |

Read `responseSummary.exists` for the dataset you named, never the list.

Why (from `views_boolean.py`): for a locus query the beacon loops over every
dataset and checks the variant inside the query already narrowed to the named
`datasetIds`, so unnamed datasets come back false. The per-dataset list is only
built for locus queries; a parameterless query returns an empty list. Whether
the Beacon v2 specification expects unevaluated datasets to be listed at all
was not checked here (the specification's `includeResultsetResponses` parameter
controls which result sets are returned; this beacon was not tested with it).

## 7. Limits

- The API answers 429 when a client is throttled; the rate-limit settings for
  this VM are in `rebuild-steps-for-taine.md` (D1 and D2). Send few queries.
- Only GRCh38, chromosomes 20 to 22, biallelic SNVs are loaded.
- A "no" for a dataset means no member of that 65 carries the allele; it says
  nothing about people outside the 65.
