#!/usr/bin/env python3
"""
interpret_tfconf_all.py — Generalized attribution analysis for any pilot(s).

Usage:
  python interpret_tfconf_all.py --pilot ets1 csl egr1 --seed 1 --out results_dir
  python interpret_tfconf_all.py --pilot all --seed 1 --out results_dir
"""
import numpy as np, torch, pickle, json, sys
from pathlib import Path
from tqdm import tqdm
from torch_geometric.data import DataLoader
from deeppbs.nn.utils import loadDataset
from deeppbs.nn import processBatch
from deeppbs.nn.metrics import mae
sys.path.insert(0, '/project2/rohs_102/shewchuk/DeepPBS/run/models')
from model_interpret import Model

R = Path('/project2/rohs_102/shewchuk/TF-conformation')
STAGE4 = R / 'output' / 'stage4_npz'
STAGE6 = R / 'output' / 'stage6_train'

def load_config(name):
    cfg_path = STAGE6 / name / name / 'config.json'
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config: {cfg_path}")
    return json.load(open(cfg_path))

def get_ckpt(name):
    ckpt = STAGE6 / name / name / 'Model.best.tar'
    if not ckpt.exists():
        raise FileNotFoundError(f"Checkpoint: {ckpt}")
    return ckpt

def find_available_arms(pilot):
    """Auto-detect baseline, augmented, augmented_relax arms for a pilot."""
    arms = []
    baseline = f'baseline_{pilot}_fold0'
    if (STAGE6 / baseline).exists():
        arms.append(baseline)
    
    augmented = f'augmented_{pilot}_fold0'
    if (STAGE6 / augmented).exists():
        arms.append(augmented)
    
    # Try to find dnarelax arm (prefer seed 1, fallback to others)
    for seed in [1, 2, 3, 4, 5, 0]:
        relax = f'augmented_{pilot}_fold0_dnarelax_s{seed}'
        if (STAGE6 / relax).exists():
            arms.append(relax)
            break
    
    return arms

def load_pilot_data(pilot):
    """Load featurized dataset for a pilot."""
    npz_dir = STAGE4 / pilot / 'output'
    files = sorted(list(npz_dir.glob('*.npz')))
    if not files:
        raise FileNotFoundError(f"No .npz in {npz_dir}")
    
    cfg = load_config(f'baseline_{pilot}_fold0')
    d,_,_,df = loadDataset([str(f) for f in files],
        nc=cfg.get('nc',4), labels_key=cfg.get('labels_key','Y_pwm'),
        data_dir=str(STAGE4), cache_dataset=False,
        balance=cfg.get('balance','unmasked'), remove_mask=False,
        scale=False, scaler=None, pre_transform=None, feature_mask=None)
    return d, df

def compute_arm(arm, batches, out_dir):
    """Compute occlusion importance for one arm."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dev = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    cfg = load_config(arm)
    ckpt_path = get_ckpt(arm)
    
    m = Model(13, 14, condition=cfg.get('condition','prot_shape'),
              readout=cfg.get('readout','all'))
    state = torch.load(ckpt_path, map_location=dev)
    m.load_state_dict(state['model_state_dict'])
    m.to(dev)
    m.eval()
    
    print(f'  Baseline: {arm}')
    base, v, x = [], [], []
    for b in tqdm(batches, desc=arm, leave=False):
        bd = processBatch(dev, b)
        with torch.no_grad():
            p = torch.softmax(m(bd['batch'], None), dim=1).cpu().numpy()
            base.append(p)
            v.append(b.v_prot.cpu().numpy())
            x.append(b.x_prot.cpu().numpy())
    
    ei = batches[0].e_prot.cpu().numpy()
    ai = np.unique(ei[0,:]).astype(int)
    
    print(f'  Occlusion: {arm} ({len(ai)} atoms)')
    diffs = []
    for b in tqdm(batches, desc='occlusion', leave=False):
        d_atom = []
        for av in ai:
            bd = processBatch(dev, b)
            with torch.no_grad():
                pm = torch.softmax(m(bd['batch'], int(av)), dim=1).cpu().numpy()
            d_atom.append(mae(base[len(diffs)], pm))
        diffs.append(np.array(d_atom))
    
    np.savez(out_dir/f'{arm}_importance.npz',
        arm=arm, prot_atom_indices=ai, occlusion_mae=np.array(diffs),
        v_prot=np.concatenate(v,axis=0), x_prot=np.concatenate(x,axis=0),
        edge_index=ei, allow_pickle=True)
    print(f'  ✓ {arm}')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--pilot', nargs='+', default=['ets1'], help='Pilot name(s), or "all"')
    p.add_argument('--out', default=f'{R}/output/interpret_results_all')
    args = p.parse_args()
    
    if args.pilot == ['all']:
        args.pilot = ['csl','egr1','engrailed','err','ets1','foxa','lef1','nfat','runx','tbp']
    
    for pilot in args.pilot:
        print(f"\n=== {pilot.upper()} ===")
        try:
            arms = find_available_arms(pilot)
            if not arms:
                print(f"  No arms found for {pilot}, skipping")
                continue
            
            print(f"  Loading data: {pilot}")
            d, df = load_pilot_data(pilot)
            dl = DataLoader(d, batch_size=1, shuffle=False)
            batches = list(dl)
            print(f"  Loaded {len(batches)} complexes")
            
            out = Path(args.out) / pilot
            out.mkdir(parents=True, exist_ok=True)
            
            for arm in arms:
                compute_arm(arm, batches, out)
        
        except Exception as e:
            print(f"  ✗ {pilot}: {e}")
            import traceback
            traceback.print_exc()
