#!/bin/bash
# proc_source.sh — put the 3DNA + DeepPBS processing binaries on PATH and point
# X3DNA at the vendored x3dna install.
#
# Self-locating: resolves the toolchain relative to THIS file
# (TF-conformation/lib/), not the current working directory. Because X3DNA is
# an absolute path, 3DNA finds its parameter files regardless of cwd — which is
# why callers no longer have to run from DeepPBS/run/process/.
_PS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # stage4_preprocess/
_REPO_DIR="$(dirname "${_PS_DIR}")"                        # TF-conformation/
export PATH="${_REPO_DIR}/lib/dependencies/bin:${PATH}"
export X3DNA="${_REPO_DIR}/lib/x3dna-v2.3-linux-64bit/x3dna-v2.3"
