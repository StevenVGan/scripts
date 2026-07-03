# scripts/pipeline/templates/

Canonical, version-controlled templates that `scripts/setup/new_project.sh`
stamps into a new project (except `README.analysis.md`, copied by hand per §11).
Editing a template here changes what every future
project is scaffolded with — this is the single source of truth for project
boilerplate (CONVENTIONS.md §3/§4 reference these files rather than inlining
their own, now-drift-prone copies).

| File | Emitted as | Notes |
|---|---|---|
| `gitignore` | `.gitignore` | one-pattern-per-line (the multi-pattern bug is designed out) |
| `samples.tsv` | `samples.tsv` | 8-col header row; human-authored sample sheet |
| `peakcall_groups.tsv` | `peakcall_groups.tsv` | headerless; `ip_bam/control_bam/name/type` (TF\|Histone) |
| `tss_groups.tsv` | `tss_groups.tsv` | csRNA only; `ip_bam/control_bam/name` for findcsRNATSS |
| `link_sample.tsv` | `link_sample.tsv` | headerless `raw_prefix -> new_sample_name` map |
| `README.project.md` | `README.md` | `__TOKENS__` substituted by new_project.sh; free-text left as TODO |
| `README.analysis.md` | `analysis/<name>/README.md` | per-analysis stub (§11) |

Headerless TSVs are documented with leading `#` comments — the pipeline readers
skip blank and `#` lines (`4.1_peak_macs3.sh:51-52`).
