"""
Stage 3 v5: DNA-aware minimization that preserves metal-coordination geometry
through pairwise sidechain restraints (no metal in the simulated system).

Per structure:
  1. PDBFixer: standard repair. BEFORE removeHeterogens, record positions of
     all structural metal ions (Zn, Mg, Mn, Fe, Ca, Co, Ni, Cu) and identify
     which protein sidechain heavy atoms coordinate each metal (within
     METAL_COORD_CUTOFF Å). Then remove all heterogens.
  2. Add Hs (Modeller, forcefield=amber14SB+gbn2).
  3. Phase 0: hydrogen-only minimization with heavy atoms frozen (mass=0).
  4. Phase 1: vdW ramping (σ=0.1→1.0) with backbone restraints AND pairwise
     harmonic restraints between metal-coordinating sidechain atoms, at their
     ORIGINAL pairwise distances. This is the key trick: instead of trying to
     parameterize the metal (which AMBER ff14SB doesn't know), we restrain the
     "coordination cage" — every pair of sidechain heavy atoms that coordinate
     the same metal stays at its starting separation.
  5. Phase 2: standard final minimization.
  6. Strip Hs, write.

The metal itself is NOT in the OpenMM system. We use it only to identify
coordination clusters in the input structure.
"""
import os, sys, argparse, time
import numpy as np
from pdbfixer import PDBFixer
import openmm
import openmm.app as app
import openmm.unit as unit

parser = argparse.ArgumentParser()
parser.add_argument("--input-pdb", required=True)
parser.add_argument("--output-pdb", required=True)
parser.add_argument("--scratch-dir", default="/scratch1/shewchuk/deeppbs_min_tmp")
parser.add_argument("--restraint-k", type=float, default=10.0,
                    help="Backbone restraint constant (kcal/mol/Å²)")
parser.add_argument("--metal-cage-k", type=float, default=20.0,
                    help="Sidechain-pair (coordination-cage) restraint k "
                         "(kcal/mol/Å²)")
parser.add_argument("--metal-coord-cutoff", type=float, default=3.0,
                    help="Distance (Å) within which a sidechain heavy atom is "
                         "considered metal-coordinating")
parser.add_argument("--ramp-stages", type=str, default="0.1,0.3,0.5,0.7,1.0")
parser.add_argument("--steps-per-stage", type=int, default=500)
parser.add_argument("--final-iterations", type=int, default=10000)
parser.add_argument("--ignore-metals", action="store_true",
                    help="Skip metal detection and coordination-cage restraints. "
                         "Reproduces the pre-metal-patch behavior — useful for A/B "
                         "testing whether the cage restraints affect downstream "
                         "training. Even if metals are present in the input PDB, "
                         "they're treated as ordinary heterogens and removed.")
args = parser.parse_args()
ramp_stages = [float(s) for s in args.ramp_stages.split(",")]

STRUCTURAL_METALS = {"ZN", "MG", "MN", "FE", "CA", "CO", "NI", "CU"}
LIGAND_RESNAMES = {"CYS", "HIS", "ASP", "GLU", "SER", "THR", "TYR",
                   "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
LIGAND_ATOM_NAMES = {"SG", "ND1", "NE2", "OD1", "OD2", "OE1", "OE2",
                     "OG", "OG1", "OH"}

t_start = time.time()
basename = os.path.basename(args.input_pdb).replace(".pdb", "")
os.makedirs(args.scratch_dir, exist_ok=True)
os.makedirs(os.path.dirname(args.output_pdb), exist_ok=True)

# ---------- Step 1: PDBFixer + identify metal coordination ----------
fixer = PDBFixer(filename=args.input_pdb)
fixer.findMissingResidues(); fixer.missingResidues = {}
fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
fixer.findMissingAtoms(); fixer.addMissingAtoms()

positions_nm_fixer = np.array([[p.x, p.y, p.z] for p in
                               fixer.positions.value_in_unit(unit.nanometer)])
cutoff_nm = args.metal_coord_cutoff / 10.0

# Find metals and their coordinating sidechain heavy atoms BEFORE
# removeHeterogens wipes them. We store ligand atoms by (chain.id, residue.id,
# atom.name) so we can re-locate them after PDBFixer/Modeller rewrites indices.
#
# In --ignore-metals mode, we skip this entirely and let removeHeterogens
# strip metals like any other heteroatom (reproducing pre-patch behavior).
metal_positions = []
coordination_clusters = []  # list of (resname, [(chain.id, residue.id, atom.name), ...])

if args.ignore_metals:
    print(f"[{basename}] --ignore-metals set: skipping metal detection (legacy mode)")
else:
    for atom in fixer.topology.atoms():
        if atom.residue.name.strip().upper() in STRUCTURAL_METALS:
            metal_positions.append((atom.residue.name.strip().upper(),
                                    positions_nm_fixer[atom.index]))

    for m_resname, m_xyz in metal_positions:
        cluster = []
        for atom in fixer.topology.atoms():
            if atom.residue.name not in LIGAND_RESNAMES:
                continue
            if atom.name not in LIGAND_ATOM_NAMES:
                continue
            if atom.element is None or atom.element.symbol == "H":
                continue
            d = np.linalg.norm(positions_nm_fixer[atom.index] - m_xyz)
            if d < cutoff_nm:
                cluster.append((atom.residue.chain.id, atom.residue.id, atom.name))
        if len(cluster) >= 2:
            coordination_clusters.append((m_resname, cluster))

    n_metals = len(metal_positions)
    n_clusters_kept = len(coordination_clusters)
    print(f"[{basename}] Found {n_metals} structural metal ion(s); "
          f"{n_clusters_kept} coordination cluster(s) with ≥2 sidechain ligands")
    for m_resname, cluster in coordination_clusters:
        ligand_strs = [f"{c[0]}/{c[1]}/{c[2]}" for c in cluster]
        print(f"[{basename}]   {m_resname} cluster ({len(cluster)} ligands): {', '.join(ligand_strs)}")

# Now remove ALL heterogens (including metals). OpenMM won't try to
# parameterize them; we'll preserve coordination via pairwise restraints.
fixer.removeHeterogens(keepWater=False)

# ---------- Step 2: Add hydrogens ----------
forcefield = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
modeller = app.Modeller(fixer.topology, fixer.positions)
modeller.addHydrogens(forcefield)
print(f"[{basename}] After H addition: {modeller.topology.getNumAtoms()} atoms")

# Build lookup: (chain.id, residue.id, atom.name) -> current modeller atom index
modeller_idx_by_id = {}
for atom in modeller.topology.atoms():
    key = (atom.residue.chain.id, atom.residue.id, atom.name)
    modeller_idx_by_id[key] = atom.index

# Translate coordination clusters to modeller atom indices
modeller_clusters = []
for m_resname, cluster in coordination_clusters:
    indices = []
    for key in cluster:
        idx = modeller_idx_by_id.get(key)
        if idx is None:
            print(f"[{basename}] WARNING: ligand {key} not found in modeller; "
                  f"skipping this ligand atom")
            continue
        indices.append(idx)
    if len(indices) >= 2:
        modeller_clusters.append((m_resname, indices))

# Compute pairwise initial distances within each cluster
positions_nm_modeller = np.array([[p.x, p.y, p.z] for p in
                                  modeller.positions.value_in_unit(unit.nanometer)])
cage_pairs = []  # list of (idx_a, idx_b, d_eq_nm)
for m_resname, indices in modeller_clusters:
    for i in range(len(indices)):
        for j in range(i+1, len(indices)):
            d_nm = np.linalg.norm(positions_nm_modeller[indices[i]] -
                                   positions_nm_modeller[indices[j]])
            cage_pairs.append((indices[i], indices[j], d_nm))
print(f"[{basename}] Coordination-cage pairwise restraints: {len(cage_pairs)}")

# ---------- Helper: clash count ----------
def count_clashes(positions, topology, threshold_nm=0.24):
    pos = np.array([[p.x, p.y, p.z] for p in positions.value_in_unit(unit.nanometer)])
    prot, dna = [], []
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol == "H": continue
        if atom.residue.name in ("DA","DG","DC","DT"): dna.append(atom.index)
        elif atom.residue.name in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
                                    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
                                    "TYR","VAL","HID","HIE","HIP"):
            prot.append(atom.index)
    if not prot or not dna: return 0, float("nan")
    p_xyz = pos[prot]; d_xyz = pos[dna]
    dists = np.sqrt(((p_xyz[:,None,:]-d_xyz[None,:,:])**2).sum(-1))
    return int((dists < threshold_nm).sum()), float(dists.min()*10)

initial_c, initial_d = count_clashes(modeller.positions, modeller.topology)
print(f"[{basename}] Initial heavy-atom clashes: {initial_c}, min_dist: {initial_d:.2f} Å")

# ============================================================
# PHASE 0: H-only minimization (freeze heavy atoms)
# ============================================================
print(f"[{basename}] Phase 0: H-only minimization (heavies frozen)")
system_h = forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.NoCutoff,
    constraints=None, rigidWater=False,
)
for i, atom in enumerate(modeller.topology.atoms()):
    if atom.element is not None and atom.element.symbol != "H":
        system_h.setParticleMass(i, 0.0)  # freeze

integrator_h = openmm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 1.0*unit.femtosecond)
try:
    platform = openmm.Platform.getPlatformByName("CUDA")
    properties = {"CudaPrecision": "mixed"}
except Exception:
    platform = openmm.Platform.getPlatformByName("CPU")
    properties = {}
print(f"[{basename}] Platform: {platform.getName()}")

sim_h = app.Simulation(modeller.topology, system_h, integrator_h, platform, properties)
sim_h.context.setPositions(modeller.positions)
e_before = sim_h.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"[{basename}]   Pre-H-min energy: {e_before:.3e} kJ/mol")
try:
    sim_h.minimizeEnergy(maxIterations=2000)
except Exception as ex:
    print(f"[{basename}] H-only min raised: {ex}")
state_h = sim_h.context.getState(getPositions=True, getEnergy=True)
e_after = state_h.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"[{basename}]   Post-H-min energy: {e_after:.3e} kJ/mol")
relaxed_positions = state_h.getPositions()

# ============================================================
# PHASE 1: vdW-ramped min with backbone + coordination-cage restraints
# ============================================================
print(f"[{basename}] Phase 1: vdW ramping with backbone + cage restraints")
system = forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds, rigidWater=False,
)

# --- Backbone restraints ---
backbone_restraint = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
backbone_restraint.addGlobalParameter("k", args.restraint_k * unit.kilocalories_per_mole / unit.angstrom**2)
backbone_restraint.addPerParticleParameter("x0")
backbone_restraint.addPerParticleParameter("y0")
backbone_restraint.addPerParticleParameter("z0")
n_restrained = 0
for atom in modeller.topology.atoms():
    is_protein_bb = (atom.residue.name in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
                                            "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
                                            "THR","TRP","TYR","VAL","HID","HIE","HIP")
                     and atom.name in ("N","CA","C"))
    is_dna_bb = (atom.residue.name in ("DA","DG","DC","DT","DA5","DG5","DC5","DT5",
                                        "DA3","DG3","DC3","DT3")
                 and atom.name in ("P","C1'"))
    if is_protein_bb or is_dna_bb:
        p = relaxed_positions[atom.index].value_in_unit(unit.nanometer)
        backbone_restraint.addParticle(atom.index, [p[0], p[1], p[2]])
        n_restrained += 1
system.addForce(backbone_restraint)
print(f"[{basename}]   Restrained {n_restrained} backbone atoms")

# --- Coordination-cage pairwise distance restraints ---
if cage_pairs:
    cage_force = openmm.HarmonicBondForce()
    k_cage = args.metal_cage_k * unit.kilocalories_per_mole / unit.angstrom**2
    for idx_a, idx_b, d_eq_nm in cage_pairs:
        cage_force.addBond(idx_a, idx_b, d_eq_nm * unit.nanometer, k_cage)
    system.addForce(cage_force)
    print(f"[{basename}]   Added {len(cage_pairs)} cage harmonic restraints")

# --- Find NonbondedForce, save sigmas ---
nb_force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
original_sigmas = [nb_force.getParticleParameters(i)[1] for i in range(nb_force.getNumParticles())]

integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 1.0*unit.femtosecond)
sim = app.Simulation(modeller.topology, system, integrator, platform, properties)
sim.context.setPositions(relaxed_positions)

def set_sigma_scale(scale):
    for i, orig in enumerate(original_sigmas):
        charge, _, epsilon = nb_force.getParticleParameters(i)
        nb_force.setParticleParameters(i, charge, orig * scale, epsilon)
    nb_force.updateParametersInContext(sim.context)

for stage_i, scale in enumerate(ramp_stages, 1):
    set_sigma_scale(scale)
    try:
        sim.minimizeEnergy(maxIterations=args.steps_per_stage)
        st = sim.context.getState(getPositions=True, getEnergy=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        c, d = count_clashes(st.getPositions(), modeller.topology)
        print(f"[{basename}]   Stage {stage_i}/{len(ramp_stages)}: σ={scale}, PE={e:.2e}, clashes={c}, min_d={d:.2f} Å")
    except Exception as ex:
        print(f"[{basename}] Stage {stage_i} FAILED: {ex}")
        sys.exit(2)

# ============================================================
# PHASE 2: Final minimization at full vdW
# ============================================================
print(f"[{basename}] Phase 2: final minimization (max {args.final_iterations} iters)")
try:
    sim.minimizeEnergy(maxIterations=args.final_iterations)
    final_state = sim.context.getState(getPositions=True, getEnergy=True)
    final_e = final_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    final_c, final_d = count_clashes(final_state.getPositions(), modeller.topology)
    print(f"[{basename}]   Final: PE={final_e:.2e}, clashes={final_c}, min_d={final_d:.2f} Å")
except Exception as ex:
    print(f"[{basename}] Final min failed: {ex}")
    sys.exit(3)

# ============================================================
# Report coordination cage drift
# ============================================================
if cage_pairs:
    final_pos = np.array([[p.x, p.y, p.z] for p in
                          final_state.getPositions().value_in_unit(unit.nanometer)])
    drifts = []
    for idx_a, idx_b, d_eq in cage_pairs:
        d_final = np.linalg.norm(final_pos[idx_a] - final_pos[idx_b])
        drifts.append((d_final - d_eq) * 10)  # Å
    print(f"[{basename}]   Cage drift (Å): "
          f"mean={np.mean(np.abs(drifts)):.2f}, max={np.max(np.abs(drifts)):.2f}")

# ---------- Strip Hs, write ----------
out_m = app.Modeller(modeller.topology, final_state.getPositions())
out_m.delete([a for a in out_m.topology.atoms() if a.element is not None and a.element.symbol == "H"])
with open(args.output_pdb, "w") as f:
    app.PDBFile.writeFile(out_m.topology, out_m.positions, f, keepIds=True)

elapsed = time.time() - t_start
print(f"[{basename}] DONE in {elapsed:.1f}s -> {args.output_pdb}")
print(f"[{basename}] Trajectory: clashes {initial_c} -> {final_c} | min_d {initial_d:.2f} -> {final_d:.2f} Å")