#!/usr/bin/env python3
"""igm_manifest.py — build paste-ready TSVs for the UCSD IGM Illumina submission form.

The IGM form has two parts, and this tool writes one TSV for each, with no comment lines,
so the values paste straight into the core's Excel:

  <prefix>_info.tsv     the submission-info block: IGM's labels VERBATIM, in form order, one
                        "label<TAB>value" row each (paste the value column next to the labels)
  <prefix>_samples.tsv  the sample table: IGM's 12 columns verbatim, one row per library

Inputs
  --info FILE ...       one or more "field<TAB>value" TSVs, applied in order (later wins).
                        Typical: a private defaults file (submitter names/emails/phone, centre
                        memberships) then the per-submission file (date, reads, platform ...).
                        Fields are matched to IGM's labels verbatim; unknown fields are an error.
  --samples FILE        a TSV whose header contains IGM's 12 column names (extra columns ignored)
  --out-prefix PREFIX   output prefix (default: IGM_manifest)

Validation (errors stop the write): unknown info fields; Run Configuration not offered for the
Platform; sample names not unique; within a pool, i7 or i5 not unique; i5 WORKFLOW B not the
reverse complement of WORKFLOW A; non-numeric conc / volume; a Library Size that is neither a
number nor a "low-high" range; empty required cells.
Warnings (written, but printed): info fields left blank.

Stdlib only; runs on the node's system python3 (3.5). Form structure as of 2026-09; labels are
copied from the form verbatim, including its typos.
"""
import argparse, csv, re, sys
from pathlib import Path

INFO_LABELS = [  # form order; "" = blank spacer row in the form
    "Date of Sample Submission", "Institute/Company Name", "PI Name", "PI Email", "Contact Name",
    "Contact Email", "Member of Moores Cancer Center", "Member of Diabetes Research Center",
    "GI Center Project", "Contact Phone Number", "Project Number", "Task Number",
    "Funding Source Number (if sposored research)", "", "",
    "Run Bioanalyzer/Tape Station", "Perform qPCR", "PhiX %", "Platform",
    "Run Configuration (Check Instructions Tab for Compatibilty)",
    "Custom Primer? (Provide more info in comments box)",
    "Number of Lanes OR Total Reads for POOL (not per sample)",
]
REQUIRED_INFO = {"Date of Sample Submission", "Institute/Company Name", "PI Name", "PI Email", "Contact Name",
                 "Contact Email", "Platform", "Run Configuration (Check Instructions Tab for Compatibilty)",
                 "Number of Lanes OR Total Reads for POOL (not per sample)"}
# Run Type -> compatible platforms (the form's "Compatible Sequencers" column)
RUN_TYPES = {
    "10X GEX (28X10X10X 90)": {"NovaSeq X Plus", "MiSeq i100"},
    "10X ATAC (50X10X24X50)": {"NovaSeq X Plus", "MiSeq i100"},
    "SR100": {"MiSeq i100"},
    "PE50": {"NovaSeq X Plus", "MiSeq i100"},
    "PE100": {"NovaSeq X Plus"},
    "PE150": {"NovaSeq X Plus", "MiSeq i100"},
    "PE300": {"MiSeq i100"},
    "PE500": {"MiSeq i100"},
    "Custom (Provide Details in Comments)": {"NovaSeq X Plus", "MiSeq i100"},
}
SAMPLE_COLS = ["Sample Name", "Pool Name", "Library Size (bp)", "Library Prep Method", "Index 1 i7 (Name)",
               "Index 1 i7 (Sequence)", "Index 2 i5 (Name)", "Index 2 i5 (Sequence-WORKFLOW A)",
               "Index 2 i5 (Sequence-WORKFLOW B)", "Conc (ng/uL)", "Volume (ul)", "Quantification Method"]
NUMERIC = {"Conc (ng/uL)", "Volume (ul)"}
SIZE_RE = re.compile(r"^\d+(\.\d+)?(\s*-\s*\d+(\.\d+)?)?$")  # Library Size: a number or a "low-high" range
RUN_LABEL = "Run Configuration (Check Instructions Tab for Compatibilty)"
READS_LABEL = "Number of Lanes OR Total Reads for POOL (not per sample)"


def read_info(paths):
    info, errs = {}, []
    for p in paths:
        with open(str(p)) as fh:
            for ln, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    errs.append("{}:{}: expected 'field<TAB>value'".format(p, ln)); continue
                k, v = parts[0].strip(), parts[1].strip()
                if k == "" or k not in INFO_LABELS:
                    errs.append("{}:{}: unknown field {!r}".format(p, ln, k)); continue
                info[k] = v
    return info, errs


def rc(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def read_samples(path):
    with open(str(path)) as fh:
        rows = [r for r in csv.reader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    hdr = rows[0]
    missing = [c for c in SAMPLE_COLS if c not in hdr]
    if missing:
        return None, ["samples: missing column(s) {}".format(missing)]
    recs = [dict(zip(hdr, r)) for r in rows[1:] if any(x.strip() for x in r)]
    errs = []
    names = [r["Sample Name"] for r in recs]
    if len(set(names)) != len(names):
        errs.append("samples: Sample Name not unique")
    for i, r in enumerate(recs, 1):
        for c in SAMPLE_COLS:
            if not r[c].strip():
                errs.append("samples row {} ({}): empty {!r}".format(i, r["Sample Name"], c))
        for c in NUMERIC:
            try:
                float(r[c])
            except ValueError:
                errs.append("samples row {}: {} = {!r} is not numeric".format(i, c, r[c]))
        if not SIZE_RE.match(r["Library Size (bp)"].strip()):
            errs.append("samples row {}: Library Size (bp) = {!r} is not a number or a low-high range".format(i, r["Library Size (bp)"]))
        a = r["Index 2 i5 (Sequence-WORKFLOW A)"].upper()
        b = r["Index 2 i5 (Sequence-WORKFLOW B)"].upper()
        if a and b and rc(a) != b:
            errs.append("samples row {}: i5 WORKFLOW B {} is not revcomp(WORKFLOW A {})".format(i, b, a))
    pools = {}
    for r in recs:
        pools.setdefault(r["Pool Name"], []).append(r)
    for pool, rs in pools.items():
        i7 = [r["Index 1 i7 (Sequence)"].upper() for r in rs]
        i5 = [r["Index 2 i5 (Sequence-WORKFLOW A)"].upper() for r in rs]
        if len(set(i7)) != len(i7):
            errs.append("pool {!r}: i7 sequences not unique".format(pool))
        if len(set(i5)) != len(i5):
            errs.append("pool {!r}: i5 sequences not unique".format(pool))
    return recs, errs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--info", action="append", required=True, metavar="FILE")
    ap.add_argument("--samples", required=True, metavar="FILE")
    ap.add_argument("--out-prefix", default="IGM_manifest")
    a = ap.parse_args()
    info, errs = read_info(a.info)
    recs, e2 = read_samples(a.samples)
    errs += e2
    plat = info.get("Platform", "")
    run = info.get(RUN_LABEL, "")
    if run and run not in RUN_TYPES:
        errs.append("Run Configuration {!r} is not one of {}".format(run, list(RUN_TYPES)))
    elif run and plat and not any(p in plat for p in RUN_TYPES[run]):
        errs.append("Run Configuration {!r} is not offered on Platform {!r} (form lists {})".format(run, plat, sorted(RUN_TYPES[run])))
    blank = [k for k in INFO_LABELS if k and not info.get(k)]
    for k in blank:
        if k in REQUIRED_INFO:
            errs.append("required info field blank: {!r}".format(k))
    if errs:
        print("ERRORS — nothing written:")
        for e in errs:
            print("  - " + e)
        sys.exit(1)
    out_info = Path(a.out_prefix + "_info.tsv")
    out_s = Path(a.out_prefix + "_samples.tsv")
    with open(str(out_info), "w") as fh:
        for k in INFO_LABELS:
            fh.write(("\t" if k == "" else "{}\t{}".format(k, info.get(k, ""))) + "\n")
    with open(str(out_s), "w") as fh:
        fh.write("\t".join(SAMPLE_COLS) + "\n")
        for r in recs:
            fh.write("\t".join(r[c] for c in SAMPLE_COLS) + "\n")
    pools = sorted({r["Pool Name"] for r in recs})
    print("wrote {}  ({} fields; blank: {})".format(out_info, sum(1 for k in INFO_LABELS if k), ", ".join(blank) if blank else "none"))
    print("wrote {}  ({} samples in {} pool(s): {}; platform {!r}, run {!r}, reads/lanes {!r})".format(
        out_s, len(recs), len(pools), ", ".join(pools), plat, run, info.get(READS_LABEL, "")))


if __name__ == "__main__":
    main()
