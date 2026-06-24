#!/bin/bash
# process_and_predict.sh — convenience wrapper: process co-crystal PDBs into npz
# features (stage4), then run DeepPBS prediction on them (stage7).
#
# Spans two stages, so it self-locates both scripts. Run it from a working dir
# that contains input.txt, process_config.json, and an output/ dir; predictions
# are written under output/.
STAGE4_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE7_DIR="$(cd "${STAGE4_DIR}/../stage7_eval" && pwd)"

# shellcheck disable=SC1091
source "${STAGE4_DIR}/proc_source.sh"
python "${STAGE4_DIR}/process_co_crystal.py" input.txt process_config.json --no_pwm
rm -f ./*.pdb ./*.par ./*.pqr ./*.r3d ./*.dat ./*.log
cat input.txt | sed 's/pdb/npz/g' | sed 's/cif/npz/g' > predict_input.txt
python "${STAGE7_DIR}/predict.py" predict_input.txt ./output/ \
    -c "${STAGE7_DIR}/pred_configs/pred_config_deeppbs.json"
