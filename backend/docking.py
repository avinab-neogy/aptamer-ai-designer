import subprocess
from pathlib import Path

def run_docking_analysis(df, target_pdb):
    # Convert target to PDBQT
    target_pdbqt = Path(target_pdb).with_suffix(".pdbqt")
    subprocess.run([
        "prepare_receptor",
        "-r", str(target_pdb),
        "-o", str(target_pdbqt)
    ])
    
    results = []
    for seq in df["sequence"]:
        # Generate ligand PDBQT
        ligand_file = Path(f"temp/{seq}.pdbqt")
        subprocess.run([
            "rna_denovo",
            "-sequence", seq,
            "-out:file:silent", str(ligand_file)
        ])
        
        # Run docking
        subprocess.run([
            "vina",
            "--receptor", str(target_pdbqt),
            "--ligand", str(ligand_file),
            "--center_x", "0", "--center_y", "0", "--center_z", "0",
            "--size_x", "20", "--size_y", "20", "--size_z", "20"
        ])
        
        # Parse results
        with open("log.txt") as f:
            affinity = float(f.readlines()[-4].split()[0])
        
        results.append({"sequence": seq, "affinity": affinity})
    
    return pd.DataFrame(results)
