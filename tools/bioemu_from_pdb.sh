#!/bin/bash
# tools/bioemu_from_pdb.sh — extract every protein chain sequence from ONE PDB
# or mmCIF structure into <PDB>_chain<X>.fasta, then submit a SLURM GPU job that
# runs BioEmu sampling + HPacker side-chain reconstruction on each chain.
#
# Run this on the login node (it calls sbatch, and may download from RCSB); the
# GPU work happens in tools/bioemu_fasta_job.sh.
#
# Sequences come from the ATOM/HETATM records of the first MODEL (observed
# residues only, matching what stage1_bioemu's PyMOL get_fastastr extraction
# does). Common modified residues (MSE, SEP, TPO, ...) map to their parent
# amino acid; nucleic acids, waters and ligands are ignored.
#
# Accepts .pdb/.ent/.cif/.mmcif, optionally .gz-compressed; the format is
# picked from the extension (and, failing that, sniffed from the contents).
# For mmCIF the auth_* chain/residue identifiers are used when present so chain
# names match the PDB-format ones.
#
# Writes (under OUTPUT_DIR, default output/stage1_bioemu/):
#   rcsb_downloads/<PDB>.<cif|pdb>          (only when fetching from RCSB)
#   <PDB>_fasta/<PDB>_chain<X>.fasta
#   <PDB>_chain<X>_conformations/            (filled by the SLURM job)
#       topology.pdb, samples.xtc             (BioEmu backbone-only)
#       samples_sidechain_rec.pdb, .xtc       (HPacker full-atom)
#
# Usage:
#   tools/bioemu_from_pdb.sh [OPTIONS] <PDB_FILE|PDB_ID> [NUM_CONFORMATIONS] [OUTPUT_DIR]
#     NUM_CONFORMATIONS  -- per chain (default: 100)
#     OUTPUT_DIR         -- default: <repo>/output/stage1_bioemu
#   Options:
#     -F, --fetch            download the structure from RCSB; the first
#                            positional argument is then a PDB ID, not a path.
#                            A bare 4-character PDB ID that is not an existing
#                            file is fetched automatically, without this flag.
#     -f, --format cif|pdb   format to download (default: cif; the legacy pdb
#                            format does not exist for large structures).
#         --force-download   re-download even if the file is already cached.
#   Env: MIN_CHAIN_LEN (default 10) -- shorter protein chains are skipped.
#
# Examples:
#   tools/bioemu_from_pdb.sh data/folds/1wui.pdb 100
#   tools/bioemu_from_pdb.sh 4guo 100              # auto-fetches 4guo.cif
#   tools/bioemu_from_pdb.sh --fetch --format pdb 1wui 50

set -eo pipefail

FETCH=0
FORMAT=cif
FORCE_DOWNLOAD=0
POSITIONAL=()

usage() {
    cat >&2 <<'USAGE'
Usage: bioemu_from_pdb.sh [OPTIONS] <PDB_FILE|PDB_ID> [NUM_CONFORMATIONS] [OUTPUT_DIR]
  NUM_CONFORMATIONS      per chain (default: 100)
  OUTPUT_DIR             default: <repo>/output/stage1_bioemu
Options:
  -F, --fetch            download the structure from RCSB; the first positional
                         argument is then a PDB ID, not a path. A bare
                         4-character PDB ID that is not an existing file is
                         fetched automatically, without this flag.
  -f, --format cif|pdb   format to download (default: cif; the legacy pdb
                         format does not exist for large structures)
      --force-download   re-download even if the file is already cached
  -h, --help             show this message
Env: MIN_CHAIN_LEN (default 10) -- shorter protein chains are skipped.
USAGE
}

while [ $# -gt 0 ]; do
    case $1 in
        -F|--fetch)
            FETCH=1; shift ;;
        -f|--format)
            FORMAT=$2; shift 2 ;;
        --format=*)
            FORMAT=${1#*=}; shift ;;
        --force-download)
            FORCE_DOWNLOAD=1; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; POSITIONAL+=("$@"); break ;;
        -*)
            echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
        *)
            POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}"

PDB_FILE=$1
NUM_CONFS=${2:-100}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "${SCRIPT_DIR}")"
OUTPUT_DIR=${3:-"${REPO_DIR}/output/stage1_bioemu"}
MIN_CHAIN_LEN=${MIN_CHAIN_LEN:-10}

if [ -z "${PDB_FILE}" ]; then
    usage
    exit 1
fi
case "${FORMAT}" in
    cif|mmcif) FORMAT=cif ;;
    pdb|ent)   FORMAT=pdb ;;
    *) echo "ERROR: --format must be 'cif' or 'pdb', got '${FORMAT}'" >&2; exit 1 ;;
esac

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"

# A bare PDB ID (4 chars, leading digit) that is not a path we can open is
# fetched from RCSB, same as if --fetch had been given.
if [ "${FETCH}" -eq 0 ] && [ ! -f "${PDB_FILE}" ] \
   && [[ "${PDB_FILE}" =~ ^[0-9][A-Za-z0-9]{3}$ ]]; then
    FETCH=1
fi

if [ "${FETCH}" -eq 1 ]; then
    if [[ ! "${PDB_FILE}" =~ ^[0-9A-Za-z]{4}$ ]]; then
        echo "ERROR: --fetch expects a 4-character PDB ID, got '${PDB_FILE}'" >&2
        exit 1
    fi
    fetch_id=$(echo "${PDB_FILE}" | tr '[:upper:]' '[:lower:]')
    DOWNLOAD_DIR="${OUTPUT_DIR}/rcsb_downloads"
    mkdir -p "${DOWNLOAD_DIR}"
    PDB_FILE="${DOWNLOAD_DIR}/${fetch_id}.${FORMAT}"
    url="https://files.rcsb.org/download/${fetch_id}.${FORMAT}"
    if [ -s "${PDB_FILE}" ] && [ "${FORCE_DOWNLOAD}" -eq 0 ]; then
        echo "Using cached download: ${PDB_FILE}"
    else
        echo "Fetching ${url}"
        tmp="${PDB_FILE}.part"
        if command -v curl >/dev/null 2>&1; then
            curl -fsSL --retry 3 -o "${tmp}" "${url}" || { rm -f "${tmp}"; fetch_rc=1; }
        elif command -v wget >/dev/null 2>&1; then
            wget -q -O "${tmp}" "${url}" || { rm -f "${tmp}"; fetch_rc=1; }
        else
            echo "ERROR: neither curl nor wget is available to fetch ${url}" >&2
            exit 1
        fi
        if [ "${fetch_rc:-0}" -ne 0 ] || [ ! -s "${tmp}" ]; then
            rm -f "${tmp}"
            echo "ERROR: download failed: ${url}" >&2
            if [ "${FORMAT}" = "pdb" ]; then
                echo "       Large structures have no legacy PDB file; try --format cif." >&2
            fi
            exit 1
        fi
        mv "${tmp}" "${PDB_FILE}"
    fi
fi

if [ ! -f "${PDB_FILE}" ]; then
    echo "ERROR: not a file: ${PDB_FILE}" >&2
    exit 1
fi

# <PDB>.pdb / <PDB>.cif / <PDB>_whatever.cif.gz / pdb<PDB>.ent -> <PDB>
pdb_id=$(basename "${PDB_FILE}")
pdb_id=${pdb_id%%.*}
pdb_id=${pdb_id%%_*}
pdb_id=${pdb_id#pdb}

FASTA_DIR="${OUTPUT_DIR}/${pdb_id}_fasta"
mkdir -p "${FASTA_DIR}"

python3 - "${PDB_FILE}" "${pdb_id}" "${FASTA_DIR}" "${MIN_CHAIN_LEN}" <<'EOF'
import gzip
import os
import shlex
import sys

pdb_file, pdb_id, fasta_dir, min_len = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common modified residues -> parent amino acid
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "HYP": "P",
    "CSO": "C", "CSD": "C", "CME": "C", "CSS": "C", "OCS": "C",
    "MLY": "K", "M3L": "K", "KCX": "K", "ALY": "K", "LLP": "K",
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "SEC": "C", "PYL": "K", "PCA": "E", "NLE": "L",
}

opener = gzip.open if pdb_file.endswith(".gz") else open


def add_residue(chains, chain, key, aa):
    chain = chain.strip() or "0"
    residues = chains.setdefault(chain, {})
    residues.setdefault(key, aa)  # first altloc wins


def parse_pdb(f):
    """Fixed-column PDB / ent records, first MODEL only."""
    chains = {}
    for line in f:
        rec = line[:6].rstrip()
        if rec == "ENDMDL":
            break  # first model only
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = line[17:20].strip()
        aa = THREE_TO_ONE.get(resname)
        if aa is None:
            continue
        # Require a CA so ligands sharing a residue name aren't picked up
        if line[12:16].strip() != "CA":
            continue
        key = (line[22:26].strip(), line[26].strip())
        add_residue(chains, line[21], key, aa)
    return chains


def cif_tokens(line):
    """Split one mmCIF loop row, honouring ' and " quoting."""
    try:
        return shlex.split(line, comments=False, posix=True)
    except ValueError:
        return line.split()


def parse_cif(f):
    """mmCIF _atom_site loop, first pdbx_PDB_model_num only."""
    chains = {}
    cols = None       # tag name -> column index, while inside the loop
    pending = []      # tag names collected in the loop header
    in_loop = False   # inside a `loop_` header
    model = None
    for raw in f:
        line = raw.strip()
        if not line or line.startswith("#"):
            in_loop = False
            pending = []
            if cols is not None:
                break  # the _atom_site loop ended
            continue
        if line == "loop_":
            in_loop = True
            pending = []
            continue
        if line.startswith("_"):
            tag = line.split()[0]
            if in_loop and tag.startswith("_atom_site."):
                pending.append(tag.split(".", 1)[1])
            elif not (in_loop and cols is None):
                pending = []
            continue
        if pending:
            cols = {name: i for i, name in enumerate(pending)}
            pending = []
        if cols is None:
            continue

        def get(*names, default=""):
            for name in names:
                i = cols.get(name)
                if i is not None and i < len(tok):
                    v = tok[i]
                    return "" if v in (".", "?") else v
            return default

        tok = cif_tokens(line)
        if len(tok) < len(cols):
            continue  # multi-line / semicolon value: not expected in _atom_site
        if get("group_PDB", default="ATOM") not in ("ATOM", "HETATM"):
            continue
        this_model = get("pdbx_PDB_model_num")
        if model is None:
            model = this_model
        elif this_model != model:
            break  # first model only
        aa = THREE_TO_ONE.get(get("label_comp_id", "auth_comp_id").upper())
        if aa is None:
            continue
        if get("label_atom_id", "auth_atom_id").strip('"') != "CA":
            continue
        key = (get("auth_seq_id", "label_seq_id"), get("pdbx_PDB_ins_code"))
        add_residue(chains, get("auth_asym_id", "label_asym_id"), key, aa)
    return chains


def is_cif(path):
    low = path.lower()
    if low.endswith(".gz"):
        low = low[:-3]
    if low.endswith((".cif", ".mmcif")):
        return True
    if low.endswith((".pdb", ".ent")):
        return False
    with opener(path, "rt") as f:  # unknown extension: sniff the contents
        for line in f:
            if line.startswith(("data_", "loop_", "_atom_site.")):
                return True
            if line[:6].rstrip() in ("ATOM", "HETATM", "HEADER", "MODEL"):
                return False
    return False


parse = parse_cif if is_cif(pdb_file) else parse_pdb
with opener(pdb_file, "rt") as f:
    chains = parse(f)

if not chains:
    sys.exit(f"ERROR: no protein residues found in {pdb_file}")

n = 0
for chain, residues in chains.items():
    seq = "".join(residues.values())
    if len(seq) < min_len:
        print(f"  skip chain {chain}: {len(seq)} residues < MIN_CHAIN_LEN={min_len}")
        continue
    path = os.path.join(fasta_dir, f"{pdb_id}_chain{chain}.fasta")
    with open(path, "w") as out:
        out.write(f">{pdb_id}_chain{chain}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")
    print(f"  chain {chain}: {len(seq)} residues -> {path}")
    n += 1

if n == 0:
    sys.exit(f"ERROR: no protein chains >= {min_len} residues in {pdb_file}")
EOF

LOG_DIR="${REPO_DIR}/slurm_output"
mkdir -p "${LOG_DIR}"
cd "${REPO_DIR}"

echo "Submitting BioEmu+HPacker: ${pdb_id} x${NUM_CONFS} confs/chain -> ${OUTPUT_DIR}"
sbatch --job-name="bioemu_${pdb_id}" \
    --output="${LOG_DIR}/bioemu_from_pdb_%j_%x.out" \
    --error="${LOG_DIR}/bioemu_from_pdb_%j_%x.err" \
    "${SCRIPT_DIR}/bioemu_fasta_job.sh" "${FASTA_DIR}" "${NUM_CONFS}" "${OUTPUT_DIR}"
