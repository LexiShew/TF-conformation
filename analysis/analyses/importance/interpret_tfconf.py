#!/usr/bin/env python3
import numpy as np, torch, pickle, json
from pathlib import Path
from tqdm import tqdm
from torch_geometric.data import DataLoader
from deeppbs.nn.utils import loadDataset
from deeppbs.nn import processBatch
from deeppbs.nn.metrics import mae
import sys
sys.path.insert(0, '/project2/rohs_102/shewchuk/DeepPBS/run/models')
from model_interpret import Model

R = '/project2/rohs_102/shewchuk/TF-conformation'
STAGE4 = Path(R) / 'output' / 'stage4_npz'
STAGE6 = Path(R) / 'output' / 'stage6_train'

def load_config(name):
    return json.load(open(STAGE6 / name / name / 'config.json'))
def get_ckpt(name):
    return STAGE6 / name / name / 'Model.best.tar'

def load_pilot_npz(pilot):
    npz_dir = STAGE4 / pilot / 'output'
    files = sorted(list(npz_dir.glob('*.npz')))
    cfg = load_config(f'baseline_{pilot}_fold0')
    d,_,_,df = loadDataset([str(f) for f in files],
        nc=cfg.get('nc',4), labels_key=cfg.get('labels_key','Y_pwm'),
        data_dir=str(STAGE4), cache_dataset=False,
        balance=cfg.get('balance','unmasked'), remove_mask=False,
        scale=False, scaler=None, pre_transform=None, feature_mask=None)
    return d, df

def compute(arm, batches, out):
    Path(out).mkdir(parents=True,exist_ok=True)
    dev = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    cfg = load_config(arm)
    ckpt_path = get_ckpt(arm)
    
    m = Model(13, 14, condition=cfg.get('condition','prot_shape'),
              readout=cfg.get('readout','all'))
    state = torch.load(ckpt_path, map_location=dev)
    m.load_state_dict(state['model_state_dict'])
    m.to(dev)
    m.eval()
    
    print(f'Baseline: {arm}')
    base, v, x = [], [], []
    for b in tqdm(batches):
        bd = processBatch(dev, b)
        with torch.no_grad():
            p = torch.softmax(m(bd['batch'], None), dim=1).cpu().numpy()
            base.append(p)
            v.append(b.v_prot.cpu().numpy())
            x.append(b.x_prot.cpu().numpy())
    
    print(f'Occlusion: {arm}')
    ei = batches[0].e_prot.cpu().numpy()
    ai = np.unique(ei[0,:]).astype(int)
    
    diffs = []
    for bi, b in enumerate(tqdm(batches)):
        d_atom = []
        for av in ai:
            bd = processBatch(dev, b)
            with torch.no_grad():
                pm = torch.softmax(m(bd['batch'], int(av)), dim=1).cpu().numpy()
            d_atom.append(mae(base[bi], pm))
        diffs.append(np.array(d_atom))
    
    np.savez(Path(out)/f'{arm}_importance.npz',
        arm=arm, prot_atom_indices=ai, occlusion_mae=np.array(diffs),
        v_prot=np.concatenate(v,axis=0), x_prot=np.concatenate(x,axis=0),
        edge_index=ei, allow_pickle=True)
    print(f'✓ {arm}')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--pilot', default='ets1')
    p.add_argument('--out', default=f'{R}/output/interpret_results')
    p.add_argument('--arms', nargs='+',
        default=['baseline_ets1_fold0','augmented_ets1_fold0','augmented_ets1_fold0_dnarelax_s1'])
    a = p.parse_args()
    
    print(f'Loading {a.pilot}')
    d, df = load_pilot_npz(a.pilot)
    dl = DataLoader(d, batch_size=1, shuffle=False)
    batches = list(dl)
    print(f'Loaded {len(batches)} complexes')
    
    for arm in a.arms:
        compute(arm, batches, a.out)
