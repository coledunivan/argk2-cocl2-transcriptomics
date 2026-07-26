# Publishing to GitHub and Zenodo

Two deposits, in this order: GitHub first (working repository), then Zenodo
(frozen, citable archive of a specific release). The GitHub–Zenodo integration
does the second one automatically once wired up.

---

## 1. Push to GitHub

```bash
cd argk2-cocl2-transcriptomics

# One-time cleanup: remove macOS metadata, leftover empty directories from the
# reorganisation, and a scratch .git created during setup. All of these are
# already gitignored — this is cosmetic, so Finder looks as clean as the repo.
rm -rf .git
find . -name '.DS_Store' -delete
rm -f .Rhistory.bak docs/.qpcr_old_readme.bak
rm -rf data/qPCR_validation_analysis outputs/figures/Supplementary

git init
git add .
git status                      # inspect before committing — see checks below
git commit -m "Analysis code and processed data for Dunivan & Fausett (G3, 2026)"

git branch -M main
git remote add origin https://github.com/coledunivan/argk2-cocl2-transcriptomics.git
git push -u origin main
```

Create the repository on GitHub first, **without** a README, .gitignore, or license —
this repo already has all three, and letting GitHub add its own creates a merge conflict
on the first push.

### Check before committing

```bash
# No file should exceed 50 MB. Expect zero output.
find . -type f -size +50M -not -path './.git/*'

# The uncompressed TFLink CSV must be ignored, the .gz must be tracked.
git check-ignore -v data/reference/RNAseq_to_TF_Targets.csv   # -> matched by .gitignore
git ls-files data/reference/                                   # -> shows only the .gz

# No macOS cruft
git ls-files | grep -c '\.DS_Store'    # -> 0
```

If `git status` shows the 79 MB CSV as untracked-but-not-ignored, stop — the `.gitignore`
rule is not matching and pushing it will bloat the repository permanently.

### Repository settings

- **Public** — required for both Zenodo archiving and G3 review.
- Description: *Analysis code and processed data for "ARGK-2 sets a transcriptional
  stress response baseline in Caenorhabditis elegans" (G3, 2026)*
- Topics: `c-elegans`, `rna-seq`, `deseq2`, `bioinformatics`, `reproducible-research`,
  `oxidative-stress`, `gene-regulatory-network`
- GitHub will detect `CITATION.cff` and add a "Cite this repository" button automatically.

---

## 2. Wire up Zenodo

1. Sign in at [zenodo.org](https://zenodo.org) with your GitHub account.
2. Go to your profile menu → **GitHub** (or zenodo.org/account/settings/github/).
3. Find `argk2-cocl2-transcriptomics` in the repository list and flip the toggle **On**.

Do this *before* creating the release. Zenodo only archives releases published after the
toggle is enabled — it will not retroactively pick up earlier ones.

---

## 3. Cut a release

On GitHub: **Releases** → **Create a new release**.

- Tag: `v1.0.0`
- Title: `v1.0.0 — G3 submission`
- Description:

  > Analysis code and processed data as submitted to G3: Genes|Genomes|Genetics.
  > Reproduces all main and supplementary figures. Raw sequencing reads: NCBI GEO GSE333535.

Publishing the release triggers Zenodo. Within a few minutes you will have a DOI.

`.zenodo.json` in the repository root controls the deposit metadata — title, authors,
affiliations, keywords, license, and the link to GSE333535 are all pre-filled, so the
Zenodo record needs no manual editing.

---

## 4. Two DOIs — use the right one

Zenodo mints two:

- **Concept DOI** — resolves to the newest version, always. Ends in a lower number.
- **Version DOI** — points at v1.0.0 specifically, permanently.

**Cite the concept DOI in the manuscript.** If you fix a bug and cut v1.0.1 during review,
the concept DOI follows it and readers land on the corrected code; a version DOI would
freeze them on the buggy release. Both appear on the Zenodo record page; the concept DOI
is labelled "Cite all versions."

---

## 5. Update the manuscript

Substitute the concept DOI into the Data Availability statement in
`docs/DATA_AVAILABILITY.md`, then paste that paragraph into the manuscript.

The GitHub username is already filled in throughout. Only the DOI remains:

```bash
grep -rn "zenodo.XXXXXXX\|ZENODO-ID" --include="*.md" --include="*.cff" .
```

---

## If you revise during review

```bash
git add -A
git commit -m "Address reviewer comments: <what changed>"
git push
```

Then cut a new release (`v1.1.0`). Zenodo archives it automatically and the concept DOI
starts resolving to it. Nothing in the manuscript needs changing.

---

## Why not a separate Zenodo data deposit?

The whole repository is ~130 MB with the TFLink file gzipped — comfortably inside
GitHub's limits and well inside Zenodo's 50 GB per-record cap. Splitting code and data
across two deposits would mean two DOIs, two download steps for reviewers, and a
synchronisation problem every time either side changes. One archive of the complete
repository is simpler and more reproducible.

The one thing that genuinely belongs elsewhere is the raw sequencing reads, which are
already at GEO under GSE333535 — that is what GEO is for, and G3 expects it.
