#!/usr/bin/env bash
# Purge structures/deeppbs_pdbs/pdb_chains.tar from TF-conformation git history
# Run from anywhere; script cd's into the repo itself.

set -euo pipefail

REPO_DIR="$HOME/carc/lab/TF-conformation"
BACKUP_DIR="$HOME/carc/lab/TF-conformation.backup.$(date +%Y%m%d-%H%M%S)"
BIG_FILE="structures/deeppbs_pdbs/pdb_chains.tar"
REMOTE_URL="https://github.com/LexiShew/TF-conformation.git"

echo "==> Step 1/8: Sanity checks"
if [ ! -d "$REPO_DIR/.git" ]; then
    echo "ERROR: $REPO_DIR is not a git repo. Aborting."
    exit 1
fi

if ! command -v git-filter-repo >/dev/null 2>&1; then
    echo "git-filter-repo not found. Installing via pip..."
    pip install --user git-filter-repo
    # Make sure pip's user bin is on PATH for this session
    export PATH="$HOME/.local/bin:$PATH"
fi

echo "==> Step 2/8: Backing up repo to $BACKUP_DIR"
cp -r "$REPO_DIR" "$BACKUP_DIR"
echo "    Backup complete. If anything breaks, restore with:"
echo "    rm -rf $REPO_DIR && mv $BACKUP_DIR $REPO_DIR"

cd "$REPO_DIR"

echo "==> Step 3/8: Confirming the big file exists in history"
if ! git log --all --oneline -- "$BIG_FILE" | grep -q .; then
    echo "WARNING: $BIG_FILE not found in git history."
    echo "Check the exact path with:  git log --all --name-only | grep pdb_chains"
    echo "Aborting so you can verify the path."
    exit 1
fi
echo "    Found $BIG_FILE in history. Proceeding."

echo "==> Step 4/8: Ensuring .gitignore excludes the file going forward"
if ! grep -qxF "$BIG_FILE" .gitignore 2>/dev/null; then
    echo "$BIG_FILE" >> .gitignore
    git add .gitignore
    git commit -m "Ignore large pdb_chains.tar" || true
    echo "    Added to .gitignore."
else
    echo "    Already in .gitignore."
fi

echo "==> Step 5/8: Undoing any partial LFS tracking for this file"
# If a previous attempt added LFS tracking, strip it so the file stays a plain
# ignored file rather than getting re-uploaded to LFS.
if [ -f .gitattributes ] && grep -q "pdb_chains.tar" .gitattributes; then
    sed -i '/pdb_chains\.tar/d' .gitattributes
    if [ -s .gitattributes ]; then
        git add .gitattributes
    else
        git rm .gitattributes
    fi
    git commit -m "Remove LFS tracking for pdb_chains.tar" || true
    echo "    Cleaned up .gitattributes."
else
    echo "    No LFS tracking to remove."
fi

echo "==> Step 6/8: Rewriting history to remove $BIG_FILE"
git filter-repo --path "$BIG_FILE" --invert-paths --force

echo "==> Step 7/8: Re-adding remote (filter-repo removes it as a safety measure)"
git remote remove origin 2>/dev/null || true
git remote add origin "$REMOTE_URL"

echo "==> Step 8/8: Force-pushing rewritten history"
git push --force --set-upstream origin main

echo ""
echo "==> Done."
echo "    The local file still exists on disk at: $REPO_DIR/$BIG_FILE"
echo "    It is now gitignored and absent from GitHub history."
echo ""
echo "Optional: reclaim local disk space by running:"
echo "    cd $REPO_DIR"
echo "    git reflog expire --expire=now --all"
echo "    git gc --prune=now --aggressive"
echo ""
echo "Once you verify the GitHub repo looks correct, you can delete the backup:"
echo "    rm -rf $BACKUP_DIR"