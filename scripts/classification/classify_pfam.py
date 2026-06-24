import numpy as np
import pandas as pd
import requests
import time

# --- CONFIGURATION ---
INPUT_FILE = 'jaspar_h11mo_cluster_wise_dna_containing_dataset.npy'
OUTPUT_FILE = 'deeppbs_tf_pfam_metadata.csv'

# Caches to prevent redundant network calls
jaspar_cache = {}
uniprot_cache = {}

def get_metadata_from_jaspar(jaspar_id):
    """Fetches UniProt ID and TF classification from JASPAR."""
    base_id = jaspar_id.split('.')[0]
    url = f"https://jaspar.elixir.no/api/v1/matrix/{base_id}/"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            # Get first UniProt ID in the list
            uniprot_id = data.get("uniprot_ids", [""])[0]
            tf_class = ", ".join(data.get("class", ["Unknown"]))
            tf_family = ", ".join(data.get("family", ["Unknown"]))
            return uniprot_id, tf_class, tf_family
    except Exception:
        pass
    return "", "Unknown", "Unknown"

def get_pfam_from_uniprot(uniprot_id):
    """Fetches Pfam accessions and names from UniProt."""
    if not uniprot_id:
        return "Unknown", "Unknown"
    
    # UniProt entry name (TAF1_HUMAN) or Accession (P21675) both work
    url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{uniprot_id}&format=json"
    if "_" in uniprot_id: # It's an entry name
        url = f"https://rest.uniprot.org/uniprotkb/search?query=id:{uniprot_id}&format=json"
        
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            results = response.json().get('results', [])
            if not results: return "Unknown", "Unknown"
            
            xrefs = results[0].get('uniProtKBCrossReferences', [])
            pfam_ids = []
            pfam_names = []
            
            for ref in xrefs:
                if ref.get('database') == 'Pfam':
                    pfam_ids.append(ref.get('id', 'Unknown'))
                    # Pfam name is usually in the first property
                    properties = ref.get('properties', [])
                    for prop in properties:
                        if prop.get('key') == 'EntryName':
                            pfam_names.append(prop.get('value'))
            
            return "|".join(pfam_ids) or "Unknown", "|".join(pfam_names) or "Unknown"
    except Exception:
        pass
    return "Unknown", "Unknown"

# 1. Load the dataset
print(f"Loading {INPUT_FILE}...")
data = np.load(INPUT_FILE, allow_pickle=True)
rows = []

print("Processing clusters and fetching Pfam data (this may take a while)...")
for cluster_idx, cluster in enumerate(data):
    for entry in cluster:
        pdb_chain = entry[0]
        pwm_ids = entry[1]
        pdb_id, chain = pdb_chain.split('_') if '_' in pdb_chain else (pdb_chain, "Unknown")
        
        for pwm_id in pwm_ids:
            uniprot_id = ""
            tf_class, tf_family = "Unknown", "Unknown"
            
            # 2. Get UniProt ID based on source
            if "jaspar" in pwm_id.lower() or pwm_id.startswith("MA"):
                if pwm_id not in jaspar_cache:
                    jaspar_cache[pwm_id] = get_metadata_from_jaspar(pwm_id)
                uniprot_id, tf_class, tf_family = jaspar_cache[pwm_id]
            elif "H11MO" in pwm_id:
                # HOCOMOCO IDs like TAF1_HUMAN.H11MO.0.A
                uniprot_id = pwm_id.split('.')[0]
            
            # 3. Get Pfam Data
            if uniprot_id not in uniprot_cache:
                print(f"  Querying Pfam for: {uniprot_id}...")
                uniprot_cache[uniprot_id] = get_pfam_from_uniprot(uniprot_id)
                time.sleep(0.1) # Respect API limits
            
            pfam_id, pfam_name = uniprot_cache[uniprot_id]
            
            rows.append({
                "PDB_ID": pdb_id.lower(),
                "Chain": chain,
                "PWM_ID": pwm_id,
                "UniProt_ID": uniprot_id,
                "TF_Class": tf_class,
                "TF_Family": tf_family,
                "Pfam_ID": pfam_id,
                "Pfam_Name": pfam_name
            })

# 4. Save to CSV
df = pd.DataFrame(rows)
df.to_csv(OUTPUT_FILE, index=False)
print(f"\nDone! Saved {len(df)} entries to {OUTPUT_FILE}")