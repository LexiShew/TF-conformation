#!/usr/bin/env bash
# finish_cleanup.sh — finish the TF-conformation repo audit items.
#
# Default: dry-run. Prints every change it would make and touches nothing.
#
#   --apply         actually perform the (reversible) file/.gitignore/.gitattributes
#                   edits and remove the superfluous local dirs
#   --migrate-lfs   ALSO run the history-rewriting `git lfs migrate import`
#                   pass. Implies --apply. Makes a timestamped backup of the
#                   repo first. You must `git push --force` afterwards.
#
# Run from anywhere; the script cd's into the repo.

set -euo pipefail

REPO_DIR="$HOME/carc/lab/TF-conformation"
APPLY=0
MIGRATE_LFS=0

while [ $# -gt 0 ]; do
    case "$1" in
        --apply)       APPLY=1 ;;
        --migrate-lfs) APPLY=1; MIGRATE_LFS=1 ;;
        -h|--help)     sed -n '2,15p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
    shift
done

[ -d "$REPO_DIR/.git" ] || { echo "Not a git repo: $REPO_DIR" >&2; exit 1; }
cd "$REPO_DIR"

log() { printf '%s\n' "$*"; }
act() { if [ "$APPLY" = 1 ]; then "$@"; else log "  [dry-run] $*"; fi; }

mode=$([ $APPLY = 1 ] && echo APPLY || echo DRY-RUN)
[ $MIGRATE_LFS = 1 ] && mode="${mode} + MIGRATE-LFS"
log "=== Mode: $mode ==="
log "Repo: $REPO_DIR"
log

###############################################################################
# 1. .gitignore — unignore lib/, dedupe pdb_chains.tar lines, add slurm_output/
###############################################################################
log "### 1. Patching .gitignore ###"

new_gitignore=$(awk '
    # Comment out lib/ and lib64/ — we use lib/ for our shared library now.
    /^lib\/$/   { print "# " $0 "  # (was for venv-style build dir; we use lib/ for shared library)"; next }
    /^lib64\/$/ { print "# " $0; next }

    # Drop duplicate pdb_chains.tar lines; we keep only the canonical one
    # (structures/deeppbs_pdbs/pdb_chains.tar). LFS will handle live versions.
    /^deeppbs_pdbs\/pdb_chains\.tar$/ { next }

    { print }

    END {
        # Add slurm_output ignore if absent.
        # (handled below via grep -q instead, since awk cant easily check
        # the whole stream for prior presence within END here.)
    }
' .gitignore)

if ! grep -qxF 'slurm_output/' .gitignore; then
    new_gitignore="$new_gitignore
# Regenerable SLURM logs
slurm_output/"
fi

if [ "$APPLY" = 1 ]; then
    cp -a .gitignore .gitignore.bak
    printf '%s\n' "$new_gitignore" > .gitignore
    log "  Patched .gitignore (backup: .gitignore.bak)"
else
    log "  [dry-run] would patch .gitignore (lib/ commented, duplicates dropped, slurm_output/ added)"
fi

###############################################################################
# 2. .gitattributes — replace 102 per-path rules with glob patterns
###############################################################################
log
log "### 2. Rewriting .gitattributes with glob patterns ###"

read -r -d '' NEW_GITATTRIBUTES <<'EOF' || true
# LFS patterns — store all large binary/structural blobs out of the regular pack.
# Code, configs, and tabular data (CSV/TSV/JSON) stay in normal git so diffs work.

# Protein/DNA structures
*.pdb     filter=lfs diff=lfs merge=lfs -text
*.cif     filter=lfs diff=lfs merge=lfs -text
*.cif.gz  filter=lfs diff=lfs merge=lfs -text
*.pqr     filter=lfs diff=lfs merge=lfs -text
*.psw     filter=lfs diff=lfs merge=lfs -text
*.pse     filter=lfs diff=lfs merge=lfs -text

# Trajectories
*.xtc     filter=lfs diff=lfs merge=lfs -text
*.dcd     filter=lfs diff=lfs merge=lfs -text

# Numerical / serialized arrays
*.npy     filter=lfs diff=lfs merge=lfs -text
*.npz     filter=lfs diff=lfs merge=lfs -text
*.pickle  filter=lfs diff=lfs merge=lfs -text
*.pkl     filter=lfs diff=lfs merge=lfs -text

# Bundles
*.tar     filter=lfs diff=lfs merge=lfs -text
*.tar.gz  filter=lfs diff=lfs merge=lfs -text
*.zip     filter=lfs diff=lfs merge=lfs -text

# Reports / large docs
*.pdf     filter=lfs diff=lfs merge=lfs -text
EOF

if [ "$APPLY" = 1 ]; then
    if [ ! -f .gitattributes.path-rules.bak ]; then
        cp -a .gitattributes .gitattributes.path-rules.bak
        log "  Backed up old .gitattributes -> .gitattributes.path-rules.bak"
    else
        log "  Existing backup .gitattributes.path-rules.bak kept (not overwritten)"
    fi
    printf '%s\n' "$NEW_GITATTRIBUTES" > .gitattributes
    log "  Wrote pattern-based .gitattributes"
else
    log "  [dry-run] would back up old .gitattributes and write pattern-based version"
    log "  [dry-run] new content preview:"
    printf '%s\n' "$NEW_GITATTRIBUTES" | sed 's/^/      /' | head -10
    log "      ... (and the rest)"
fi

###############################################################################
# 3. Remove superfluous local dirs (and drop from index if tracked)
###############################################################################
log
log "### 3. Removing superfluous local paths ###"

declare -a SUPERFLUOUS=(
    "structures/deeppbs_pdbs/little_test"
    "scripts/__pycache__"
    "slurm_output"
)

for p in "${SUPERFLUOUS[@]}"; do
    if [ -e "$p" ]; then
        if git ls-files --error-unmatch "$p" >/dev/null 2>&1 \
           || git ls-files --error-unmatch "$p/" >/dev/null 2>&1 ; then
            log "  git rm -rf --cached '$p'  (tracked — keep working tree, drop from index)"
            act git rm -rf --cached "$p" -q 2>/dev/null || true
        fi
        log "  rm -rf '$p'"
        act rm -rf "$p"
    else
        log "  [gone already] $p"
    fi
done

###############################################################################
# 4. Stage lib/ now that it's no longer ignored
###############################################################################
log
log "### 4. Staging lib/ (no longer ignored) ###"

if [ "$APPLY" = 1 ]; then
    # Force re-evaluation of cached ignores; safe even if nothing changes.
    log "  git add .gitignore .gitattributes"
    act git add .gitignore .gitattributes
    log "  git add lib/"
    act git add lib/
    log "  (commit yourself with a message like: 'cleanup: unignore lib/, glob LFS rules, drop junk')"
else
    log "  [dry-run] would: git add .gitignore .gitattributes lib/"
fi

###############################################################################
# 5. (Optional) LFS history migration
###############################################################################
if [ $MIGRATE_LFS = 1 ]; then
    log
    log "### 5. LFS history migration (DESTRUCTIVE — rewrites history) ###"

    BACKUP_DIR="$HOME/carc/lab/TF-conformation.backup.$(date +%Y%m%d-%H%M%S)"
    log "  Backing up to: $BACKUP_DIR"
    cp -a "$REPO_DIR" "$BACKUP_DIR"

    if ! command -v git-lfs >/dev/null 2>&1; then
        echo "git-lfs not installed; aborting migration." >&2
        exit 1
    fi
    if ! command -v git-filter-repo >/dev/null 2>&1; then
        log "  (git-filter-repo not strictly required for 'lfs migrate', but recommended)"
    fi

    log "  Running: git lfs migrate import --include='*.pdb,*.cif,*.cif.gz,*.pqr,*.pse,*.psw,*.xtc,*.dcd,*.npy,*.npz,*.pickle,*.pkl,*.tar,*.tar.gz,*.zip,*.pdf' --everything"
    git lfs migrate import \
        --include='*.pdb,*.cif,*.cif.gz,*.pqr,*.pse,*.psw,*.xtc,*.dcd,*.npy,*.npz,*.pickle,*.pkl,*.tar,*.tar.gz,*.zip,*.pdf' \
        --everything

    log
    log "  Migration done. Verify with:  git lfs ls-files | head"
    log "  Then push (force, since history was rewritten):  git push --force --all && git push --force --tags"
    log "  Restore from backup if anything looks wrong:"
    log "    rm -rf '$REPO_DIR' && mv '$BACKUP_DIR' '$REPO_DIR'"
fi

log
log "Done."
[ $APPLY = 1 ] || log "(dry-run; re-run with --apply to act, or --migrate-lfs to also rewrite history)"
