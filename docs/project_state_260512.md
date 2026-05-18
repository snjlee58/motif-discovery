# Project state — advisor meeting prep (2026-05-12)

## TL;DR (2 sentences)

Pipeline produces bimodal F1 on M-CSA monomers — F1=1 on ~9 PDBs, F1=0 on ~108 PDBs, with the rest spread between. After a week of diagnostics, the root cause is that the column-level conservation signal (AA *and* 3Di) **does not differentiate catalytic residues for the failure cohort**, and the most promising rescue I could think of (lift-based scoring) made things worse, not better.

## What the pipeline does

| Step | Script | What it does |
|---|---|---|
| 0 | `pipeline.sh` | Resolve PDB → UniProt via UniProt ID Mapping API |
| 1 | `pipeline.sh` | Find AFDB cluster representative; pull all AFDB50 reps (`cluFlag=1`) — this is the ~50%-identity-deduped homolog list |
| 2 | `pipeline.sh` | Parallel `wget` AlphaFold PDBs for all cluster members + RCSB experimental PDB for the query |
| 3 | `pipeline.sh` (foldmason) | Structural MSA via FoldMason → `foldmason_result_aa.fa` and `foldmason_result_3di.fa` |
| 4 | `src/score_conservation.py` | Per-column conservation scoring (AA = 50% inverted Shannon entropy + 50% BLOSUM62 self-similarity; 3Di = inverted Shannon entropy only). Outputs per-column `conservation`, `identity`, `gap_frequency`, `3di_conservation`. |
| 5 | `src/map_alignment_to_pdb.py` | Map alignment columns to PDB `auth_resid` using the query's PDB ATOM records |
| 5b | (P2Rank) | Optional binding-pocket prediction; non-fatal if unavailable |
| 6 | `src/extract_top_conserved.py` | Pick top-N most-conserved residues (sanity-check output, not used for benchmark) |
| 7 | `src/benchmark_mcsa.py` | Score every column with the full formula (below), select top-N, compute precision/recall/F1 vs M-CSA ground truth |

## Scoring formula (the heart of the project)

```
base_score = w_cons · conservation
           + w_p2rank · P2Rank_pocket_probability
           + w_3di · 3Di_structural_conservation

Pass 2: final_score = base_score + w_cluster · (count of high-scoring neighbors within 15Å)
```

Current weights (hand-tuned, in `src/benchmark_mcsa.py`):

| Weight | Value | Notes |
|---|---|---|
| `w_cons` | 1.0 | dominant |
| `w_p2rank` | 0.35 | binding-pocket prior |
| `w_3di` | 0.30 | 3Di entropy at the column |
| `w_cluster` | 0.30 | spatial clustering bonus |

**Note on the catalytic-propensity term (removed 2026-05-12):** Earlier versions of the formula included a `w_prop = 0.25` term using per-AA catalytic frequencies from Ribeiro et al. 2017 (His = 8.01, Cys = 4.66, etc.). I removed this prior because:
1. It's computed *from M-CSA itself*, so benchmarking against M-CSA with it active is circular — the score inflates the F1 reported on M-CSA.
2. It doesn't generalise to non-classical active sites (DNA/RNA binding, ligand interactions) — it's specific to enzymatic catalysis as M-CSA defines it.

Pre-removal benchmark numbers (mean F1 = 0.313 on M-CSA monomers, the bimodal result discussed below) included the propensity term. **The benchmark needs to be re-run without it to know the true baseline.** Likely lowers the F1 somewhat for the win cohort but doesn't change the bimodality (the failure cohort gets 0 either way).

Selection: `--top-n auto` (predict N = number of M-CSA ground-truth residues — same convention Squidly uses). Filters: `--exclude-gaps` (drop columns with gap_frequency > 0.5), `--min-identity 0.2` (drop columns where consensus AA appears in <20% of sequences).

## Benchmark setup

- **PDB list**: `benchmark/mcsa_representatives_parsed_monomers.tsv` — M-CSA representatives filtered to monomeric chains (354 PDBs).
- **Ground truth**: `$FAST/m-csa/catalytic_residues_homologues_parsed.tsv` — M-CSA *homologues* file (includes inferred catalytic residues for homologs, not just experimentally validated).
- **Compute**: 302/354 PDBs completed end-to-end via SLURM batch (`scripts/motif_benchmark.sbatch` → `benchmark/run_batch.sh`).
- **Mean F1**: 0.313 on completed set; **F1 distribution is bimodal**.

## What I tested (in order)

| Diagnostic | Script | What it answered |
|---|---|---|
| Top-N curve at 1×/2×/3× ground truth | `analysis/analyze_top_n_comparison.py` | 74/277 PDBs still F1=0 at 3× → residues genuinely absent from top-3N |
| Mapping vs scoring split | `analysis/diagnose_missing_residues.py` | Of those 74, only 5 are mapping failures; 69 are scoring failures (catalytic residues ARE in `alignment_mapping.json`) |
| Filter vs scoring split | `analysis/diagnose_filter_vs_scoring.py` | 63% of scoring-failure residues are filter losses (`exclude-gaps` + `min-identity` together); 37% pass filters but rank too low |
| 3Di rescue test | (same script, 3Di rank added) | **0/332** catalytic residues are caught by 3Di rank when AA rank misses them |
| Filter relaxation rerun | `EXCLUDE_GAPS=0 MIN_IDENTITY=0 ... rescore_batch.sh` | Mean F1: 0.325 → 0.313 (slight DROP). Filters were not the bottleneck. |
| Wins-vs-fails MSA quality | `analysis/diagnose_msa_quality.py` | Wins: cat_identity 0.85, lift +0.77, MSA depth 1347 sequences. Fails: cat_identity 0.30, lift +0.23, MSA depth 3000 sequences |
| Lift-based re-ranker | `analysis/test_lift_scoring.py` | Lift_identity ranking destroys wins (1.0 → 0.39) and barely rescues fails (10/108 to F1>0). Net F1: 0.313 → 0.201 |

## What I found

### The bimodality is intrinsic to the data, not the scoring

Across every diagnostic, the failure cohort has the same pattern: catalytic residues sit at columns that are **not at the top** of any column-level conservation signal. AA conservation, AA lift, 3Di conservation, 3Di lift, combinations — none of them rank catalytic columns above non-catalytic ones in the failure PDBs.

### Wins are easy

When catalytic residues *are* sequence-conserved (cat_identity 0.85 — nearly universal), any reasonable scoring catches them. The pipeline already achieves mean F1 ≈ 0.78 on the high+perfect cohorts (~109 PDBs).

### 3Di does not add a column-level signal

The original framing was "3Di conservation catches catalytic residues that AA conservation misses." The data robustly contradicts this:
- 0/332 catalytic columns rescued by 3Di rank when AA rank missed them
- 3Di-only ranking achieves F1 ≈ 0.00 on the perfect cohort — it doesn't even catch the easy cases
- 3Di and AA conservation are tightly correlated; 3Di has smaller dynamic range across columns and effectively re-ranks the same way AA does

### MSA depth correlates *negatively* with F1

Failure-cohort MSAs have 3000 sequences on average; success-cohort MSAs have 1347. Deeper MSAs of distant homologs *dilute* the catalytic signal — many non-catalytic columns become more conserved than catalytic ones once the alignment is broad enough.

### Filter relaxation does not help

I built the hypothesis that `--exclude-gaps` + `--min-identity 0.2` were dropping catalytic columns. After plumbing them into env vars and rerunning rescore with both disabled: mean F1 stays at 0.313 (slight drop, not improvement). Filters were not the issue.

## Root cause, best current understanding

**For ~36% of M-CSA monomer test PDBs, the catalytic residues are evolutionarily flexible — they vary across the AFDB cluster more than non-catalytic structural residues vary.** No column-level conservation signal will rank them at the top.

This is biologically plausible: many catalytic residues serve a *role* (general acid/base, nucleophile, electrostatic stabilizer) that can be filled by multiple amino acids. D and E can both be acid/base; H and K can both be a base; C and S can both be a nucleophile. Conservation methods only catch the residues that are universally one specific amino acid.

## Caveats / unresolved confounders

1. **Ground truth source**: M-CSA's homologues file (`catalytic_residues_homologues_parsed.tsv`) contains residues inferred from related enzymes, not all experimentally validated. There may be label noise — some "F1=0" results could be the pipeline correct but M-CSA wrong.
2. **PDB list source vs ground truth**: PDBs come from `mcsa_representatives_parsed_monomers.tsv` (the representatives file, filtered to monomers); ground truth uses the homologues file. Schema mismatch is possible.
3. **No independent benchmark yet**: haven't tested against Squidly's CataloDB or the 6 legacy datasets. Currently downloading. Until that's done, I can't separate "pipeline is bad" from "M-CSA benchmark is noisy."

## Options I'm weighing

| Option | What I'd do | Cost | Best case | Worst case |
|---|---|---|---|---|
| **A. MSA sub-sampling** | Cap `cluster_members.txt` at top-K closest homologs (currently uses all up to ~3000). Wins have 1347 sequences; failures have 3000. Maybe shallow MSAs concentrate signal. | One pipeline.sh modification; rerun on failure-cohort PDBs (~108) | F1 lift on the failure cohort if depth dilutes signal | No change, confirms the failures are intrinsically conservation-blind |
| **B. Honest reframe** | Pivot from "novel signal" to "characterize the conservation-prediction boundary." Stratify results: pipeline achieves F1 ≈ 0.78 on the ~60% of enzymes where catalytic residues are evolutionarily conserved; characterize the failure boundary biologically. | Mostly analysis + writing | Defensible contribution even without a new method | Reviewers ask "why not just use Squidly?" |
| **C. Squidly benchmarks** | Run on CataloDB (232 sequences, <30% identity, the harder benchmark Squidly emphasizes) and the 6 legacy datasets. In progress — data is downloading. | One conversion script + a SLURM batch on ~232 new PDBs | Validates whether bimodality is M-CSA-specific or general; gives a direct SOTA comparison | If pipeline F1 is much lower than Squidly's on the same data, that's a real problem |
| **D. Add a non-conservation signal** | The failure cohort is conservation-blind. Need something else: P2Rank weight up, learned per-PDB cluster geometry, ML-supervised scoring on Squidly's training data. | Moderate to large | Could close the gap on the failure cohort | "Supervised scoring" trades away the "no training" property — may not be a clean novelty |

## What I'd want guidance on

1. **Is "characterize the boundary" defensible as a thesis contribution?** I have a method that achieves F1 ≈ 0.78 on enzymes with evolutionarily-conserved catalytic residues, F1 ≈ 0 otherwise, and a clear biological explanation for why. Is that worth writing up on its own, or do I need a fix for the failure cohort to publish?

2. **MSA sub-sampling — worth trying or noise?** The depth correlation (3000 fail vs 1347 win) is striking. Is there reason to believe shallower MSAs would have a stronger conservation signal? Or is depth a *consequence* of broad enzyme families (which fail), not a *cause*?

3. **3Di is robustly column-blind in my data.** Is there a way to use 3Di that isn't column-by-column ranking? (Structural cluster filtering? Distance-based features between conserved positions? Local 3Di patches?) I want to keep some structural-conservation angle alive since that's the lab specialty.

4. **Supervised scoring**: would learning the scoring weights from Squidly's training set violate the "no training, no reaction input" claim? Or is small-scale weight tuning (a handful of weights, not deep learning) acceptable framing-wise?

5. **Benchmark recommendation**: should I prioritize CataloDB (low-identity, the Squidly novelty case) or the 6 legacy datasets (direct comparison to AEGAN/SCREEN/CRpred)? I have compute for one of these in the near term.

## Reproducibility / where everything lives

| What | Path |
|---|---|
| Pipeline entry point | [pipeline.sh](../pipeline.sh) |
| Scoring | [src/benchmark_mcsa.py](../src/benchmark_mcsa.py) |
| Conservation lib | [src/conservation.py](../src/conservation.py) |
| Batch run results | `$FAST/motif/batch_runs/260426_211019_job527128/` (302 PDBs) |
| Original rescore | `260426_211019_job527128/rescore_260427_162231_topn1-2-3/` |
| Filter-relaxed rescore | `260426_211019_job527128/rescore_260512_004120_topn1-2-3_nogap_id0/` |
| Diagnostics + per-PDB tables | `analysis/*.tsv` |
| Diagnostic scripts | `analysis/analyze_top_n_comparison.py`, `diagnose_missing_residues.py`, `diagnose_filter_vs_scoring.py`, `diagnose_msa_quality.py`, `test_lift_scoring.py` |
