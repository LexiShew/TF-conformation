import os
import csv
import time
import requests
import re
from pathlib import Path
from Bio import SeqIO

# --- Configuration ---
EBI_API_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
POLL_INTERVAL = 15  # Seconds to wait between status checks
OUTPUT_CSV = "dna_binding_domains.csv"
EMAIL = "shewchuk@usc.edu" # EBI requires a valid email for API tracking

def submit_job(sequence):
    """Submits sequence to InterProScan API."""
    params = {
        'email': EMAIL,
        'title': 'TF_Domains',
        'goterms': 'true',
        'stype': 'p',
        'sequence': sequence,
        'appl': 'PfamA'
    }
    try:
        response = requests.post(f"{EBI_API_URL}/run", data=params, timeout=30)
        response.raise_for_status()
        return response.text
    except Exception as e:
        print(f"  [!] Error submitting job: {e}")
        return None

def wait_for_job(job_id):
    """Polls the API until the job is finished."""
    while True:
        try:
            status_res = requests.get(f"{EBI_API_URL}/status/{job_id}", timeout=20)
            status = status_res.text
            if status == "FINISHED":
                return True
            if status in ["FAILURE", "NOT_FOUND"]:
                return False
        except:
            pass
        time.sleep(POLL_INTERVAL)

def get_dna_binding_hits(job_id, original_sequence):
    """Fetches results and filters for DNA-binding Pfam domains."""
    try:
        res = requests.get(f"{EBI_API_URL}/result/{job_id}/json", timeout=30)
        res.raise_for_status()
        data = res.json()
    except:
        return []

    hits = []
    # InterProScan JSON structure: results[0] -> matches
    for match in data.get('results', [{}])[0].get('matches', []):
        signature = match.get('signature', {})
        entry = signature.get('entry', {})
        
        # Check metadata for DNA binding markers
        name = (entry.get('name') or "").lower()
        description = (entry.get('description') or "").lower()
        go_terms = [go.get('id') for go in entry.get('goGraphNodes', [])]
        
        if "dna-binding" in name or "dna binding" in description or "GO:0003677" in go_terms:
            for loc in match.get('locations', []):
                start, end = loc['start'], loc['end']
                hits.append({
                    "pfam_id": signature.get('accession'),
                    "domain_name": entry.get('name', 'Unknown'),
                    "description": entry.get('description', 'N/A'),
                    "start": start,
                    "end": end,
                    "seq_fragment": original_sequence[start-1:end]
                })
    return hits

def process_directory(base_dir):
    base_path = Path(base_dir)
    
    # Initialize CSV file
    file_exists = os.path.isfile(OUTPUT_CSV)
    with open(OUTPUT_CSV, 'a', newline='') as csvfile:
        fieldnames = ['pdb_id', 'pfam_id', 'domain_name', 'start', 'end', 'sequence', 'description']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()

        # Iterate through the structure: <DIR>/<PDB>_chains/<PDB>_conformations/sequence.fasta
        # Using glob to find all matching fasta files
        fasta_files = list(base_path.glob("*_chains/*_conformations/sequence.fasta"))
        
        print(f"Found {len(fasta_files)} sequences to process.")

        for fasta_path in fasta_files:
            # Extract PDB ID from the parent folder name (e.g., '1a1f_chains' -> '1a1f')
            pdb_id = fasta_path.parts[-3].split('_')[0]
            print(f"Processing PDB: {pdb_id}...")

            # Load sequence
            try:
                record = SeqIO.read(fasta_path, "fasta")
                seq_str = str(record.seq)
            except Exception as e:
                print(f"  [!] Failed to read {fasta_path}: {e}")
                continue

            # API Workflow
            job_id = submit_job(seq_str)
            if job_id and wait_for_job(job_id):
                domain_hits = get_dna_binding_hits(job_id, seq_str)
                
                if not domain_hits:
                    print(f"  [-] No DNA-binding domains found for {pdb_id}")
                
                for hit in domain_hits:
                    writer.writerow({
                        'pdb_id': pdb_id,
                        'pfam_id': hit['pfam_id'],
                        'domain_name': hit['domain_name'],
                        'start': hit['start'],
                        'end': hit['end'],
                        'sequence': hit['seq_fragment'],
                        'description': hit['description']
                    })
                    csvfile.flush() # Ensure it writes to disk immediately
                print(f"  [+] Success: {len(domain_hits)} domains found.")
            else:
                print(f"  [!] Job failed for {pdb_id}")

if __name__ == "__main__":
    # Replace '.' with the actual path to your <DIR>
    target_dir = input("Enter the path to the root directory (<DIR>): ").strip()
    if os.path.isdir(target_dir):
        process_directory(target_dir)
    else:
        print("Invalid directory.")