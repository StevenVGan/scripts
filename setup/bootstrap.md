# bootstrap.md — stand up a new analysis node

Self-contained runbook for the new box. Assumes internet + a shell. Companion to
the private `work_brain/PROVISION.md` (fuller rationale + Phase 2 references + §5 sync).

## 1. Base tooling
1. **miniforge** at `$HOME/miniforge3` (default location — keeps
   `CONDA_BIO_ENV=${HOME}/miniforge3/envs/bio` valid). Install `mamba`.
2. **gh** — install fresh (do NOT copy the primary box's binary; it's
   glibc-2.23-linked). Then:
   ```bash
   gh auth login --hostname github.com --git-protocol https   # account StevenVGan
   gh auth setup-git            # points git's credential helper at the new gh
   gh config set git_protocol https
   git config --global user.name 'Steven Gan'
   git config --global user.email digan@ucsd.edu
   ```

## 2. Clone the spine
```bash
# copy clone_all.sh over (or clone scripts first, then use scripts/setup/clone_all.sh)
WORK=$HOME/work bash clone_all.sh            # scripts/ (public) + work_brain/ (private brain)
#   add --all / --tf-atlas / --pwm / --idr to also clone the standalone tool repos
```

## 3. Rebuild envs (bio + sc)
```bash
mamba env create -f $HOME/work/scripts/env/bio.yml
mamba env create -f $HOME/work/scripts/env/sc.yml
```
Fresh solve, NOT conda-pack. On a newer-glibc box you should NOT need
`CONDA_OVERRIDE_GLIBC`. If a full-pin solve fails, retry with `--no-builds`, then
re-freeze. Watch `ucsc-liftover` / `bismark` (historic zlib pin conflicts).

## 4. Drop in the carry-bundle (durable Claude brain + config)
```bash
CARRY=$HOME/work/work_brain                                    # the private brain repo (cloned above)
MEMDIR=$HOME/.claude/projects/$(printf '%s' "$HOME/work" | sed 's#/#-#g')/memory
mkdir -p "$MEMDIR"; cp "$CARRY"/memory/*.md "$MEMDIR"/          # 29 durable memories + MEMORY.md
mkdir -p $HOME/work/.claude; cp "$CARRY"/settings.local.json $HOME/work/.claude/
cp "$CARRY"/CLAUDE.md.template $HOME/work/CLAUDE.md            # then edit: fill <WORK_ROOT>/<REF_ROOT>/<N>, drop the .template line
```
Consider writing 1–2 explicit `user`-type memories (there are none yet).

## 5. Verify
```bash
gh auth status                                   # logged in to StevenVGan over https
git -C $HOME/work/scripts ls-remote origin >/dev/null && echo "git auth OK"
$HOME/work/scripts/setup/new_project.sh --help   # scaffolder present
# ref-free smoke test: peak_ops.sh on two BEDs with --genome-sizes supplied
```

## 6. First project
```bash
scripts/setup/new_project.sh cnr <Assay_date_target_cell_person>
# fill link_sample.tsv / samples.tsv / peakcall_groups.tsv, then:
cd seq/cnr/<name>/script && ./link_fastq.sh && ./run_all.sh
```
Running a pipeline end-to-end needs references — that's **Phase 2** (work_brain/PROVISION.md §4):
copy the former-NAS bowtie2 indexes + FASTAs under `$HOME/work/ref` and set `REF_ROOT`.

## Ongoing
- Pull convention/tool updates: `git -C scripts pull`.
- Pull durable memory updates: `scripts/setup/sync_brain.sh pull`.
- See work_brain/PROVISION.md §5 for the two-node sync discipline.
