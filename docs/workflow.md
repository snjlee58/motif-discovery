# Workflow

Visual reference for how the motif-discovery pipeline is invoked and what
each step does. Diagrams are written in [Mermaid](https://mermaid.js.org/) —
they render natively on GitHub and in VS Code with the *Markdown Preview
Mermaid Support* extension.

## 1. Orchestration

How a run gets started, persisted, aggregated, or re-scored.

```mermaid
flowchart TD
    User((User))
    User -->|one PDB| SR[scripts/single_run.sbatch]
    User -->|all monomers| MB[scripts/motif_benchmark.sbatch]

    SR --> PS[pipeline.sh PDB_ID]

    MB --> RB[benchmark/run_batch.sh]
    TSV[(mcsa_representatives_<br/>parsed_monomers.tsv)] -->|reads PDBs| RB
    RB -->|xargs -P N| PS

    PS --> Out[(per-PDB results dir<br/>baseline_performance.json<br/>top5_motif.txt<br/>conservation.json<br/>alignment_mapping.json<br/>...)]

    Out -.aggregate.-> SR2[benchmark/summarize_results.py]
    Out -.re-score only.-> RBO[benchmark/run_benchmark_only.sh]

    SR2 --> Table[/F1 / Precision / Recall<br/>summary table/]
    RBO --> RescoreOut[(rescore_TS/<br/>fresh JSONs + summary)]
```

## 2. Inside `pipeline.sh`

The per-PDB steps, in order.

```mermaid
flowchart TD
    Start([pipeline.sh PDB_ID]) --> S0
    S0[Step 0 · Resolve UniProt<br/>UniProt ID Mapping API] --> S1
    S1[Step 1 · Find AFDB cluster<br/>grep $FAST/afdb_clusters/...tsv<br/>→ REP_ID + member list] --> S2
    S2[Step 2 · Download structures<br/>xargs parallel wget AlphaFold<br/>+ RCSB experimental PDB] --> S3
    S3[Step 3 · Structural MSA<br/>foldmason easy-msa<br/>→ foldmason_result_aa.fa<br/>→ foldmason_result_3di.fa] --> S4
    S4[Step 4 · Conservation scoring<br/>src/score_conservation.py] --> S5
    S5[Step 5 · Alignment → PDB mapping<br/>src/map_alignment_to_pdb.py] --> S5b
    S5b[Step 5b · P2Rank binding sites<br/>prank predict<br/>→ src/parse_p2rank.py] --> S6
    S6[Step 6 · Top-N conserved positions<br/>src/extract_top_conserved.py] --> S7
    S7[Step 7 · M-CSA benchmark<br/>src/benchmark_mcsa.py] --> Done([Result files])

    Lib[(src/conservation.py<br/>library)] -.imported by.-> S4

    S4 --> J4[(_conservation.json)]
    S5 --> J5[(alignment_mapping.json)]
    S5b --> J5b[(p2rank_scores.json)]
    S6 --> J6[(top5_motif.txt<br/>top10_motif.txt)]
    S7 --> J7[(baseline_performance.json)]
```
