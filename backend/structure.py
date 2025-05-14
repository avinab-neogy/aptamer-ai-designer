import os
import subprocess
import hashlib
import numpy as np
import shutil
from pathlib import Path
import RNA  # ViennaRNA package

# Configuration
TEMP_DIR = Path("temp")
TEMP_DIR.mkdir(exist_ok=True, parents=True)

# SimRNA paths
SIMRNA_PATH = Path.home()/"Downloads"/"SimRNA_64bitIntel_Linux"
SIMRNA_BIN = SIMRNA_PATH/"SimRNA"
SIMRNA_TRAFL2PDBS = SIMRNA_PATH/"SimRNA_trafl2pdbs"
SIMRNA_DATA = SIMRNA_PATH/"data"

# Create data symlink in working directory
DATA_LINK = Path("data")
if not DATA_LINK.exists():
    try:
        DATA_LINK.symlink_to(SIMRNA_DATA, target_is_directory=True)
        print("Created SimRNA data symlink")
    except Exception as e:
        print(f"Failed to create data symlink: {str(e)}")

def predict_secondary_structure(sequence):
    """Predict RNA secondary structure using ViennaRNA"""
    rna_seq = sequence.replace('T', 'U')
    structure, mfe = RNA.fold(rna_seq)
    return structure, mfe

def generate_dummy_pdb(sequence, output_path):
    """Generate dummy PDB structure"""
    with open(output_path, 'w') as f:
        f.write("HEADER    SIMULATED APTAMER STRUCTURE\n")
        f.write("TITLE     APTAMER VISUALIZATION MODEL\n")
        
        # Create a simple 3D helix structure
        x, y, z = 0, 0, 0
        helix_radius = 5.0
        rise_per_base = 3.4
        bases_per_turn = 10
        
        for i, base in enumerate(sequence):
            angle = (i / bases_per_turn) * 2 * np.pi
            x = helix_radius * np.cos(angle)
            y = helix_radius * np.sin(angle)
            z = i * rise_per_base
            
            f.write(f"ATOM  {i+1:5d}  C{base}  {base} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
            if i > 0:
                f.write(f"CONECT{i:5d}{i+1:5d}\n")
        
        f.write("END\n")
    return output_path

def predict_aptamer_structure(sequence):
    """Predict 3D structure using SimRNA with proper trajectory handling"""
    rna_sequence = sequence.replace('T', 'U').lower()
    ss, mfe = RNA.fold(rna_sequence)
    
    # Create unique working directory
    session_id = hashlib.md5(sequence.encode()).hexdigest()[:8]
    work_dir = TEMP_DIR/f"simrna_{session_id}"
    work_dir.mkdir(exist_ok=True, parents=True)  # Added parents=True
    
    # File paths
    seq_file = work_dir/"input.seq"
    ss_file = work_dir/"input.ss"
    output_prefix = work_dir/"output"
    output_pdb = TEMP_DIR/f"aptamer_{session_id}.pdb"
    
    try:
        # Write input files (validate first)
        seq_file.write_text(rna_sequence)
        ss_file.write_text(ss)
        print(f"Input files created: {seq_file.exists()}, {ss_file.exists()}")

        # 1. Generate initial reference structure with verbose output
        print("Running initial SimRNA (0 iterations)...")
        result_init = subprocess.run([
            str(SIMRNA_BIN),
            "-s", str(seq_file),
            "-S", str(ss_file),
            "-c", str(SIMRNA_PATH/"configSA.dat"),
            "-o", str(output_prefix),
            "-n", "0"
        ], capture_output=True, text=True, check=True)
        
        print(f"Initial run stdout: {result_init.stdout}")
        print(f"Initial run stderr: {result_init.stderr}")
        print(f"Files after initial run: {os.listdir(work_dir)}")

        # Find initial PDB with multiple fallbacks
        init_pdb = next(work_dir.glob("*-000001.pdb"), None)
        if not init_pdb:
            init_pdb = next(Path().glob("*-000001.pdb"), None)  # Check current dir
        if not init_pdb:
            raise FileNotFoundError(f"Initial PDB not found in {work_dir} or current directory")
        
        print(f"Found initial PDB at: {init_pdb}")

        # 2. Run full simulation with output capture
        print("Running full simulation...")
        result_sim = subprocess.run([
            str(SIMRNA_BIN),
            "-s", str(seq_file),
            "-S", str(ss_file),
            "-c", str(SIMRNA_PATH/"configSA.dat"),
            "-o", str(output_prefix),
            "-n", "10000"
        ], capture_output=True, text=True, check=True)
        
        print(f"Simulation stdout: {result_sim.stdout}")
        print(f"Simulation stderr: {result_sim.stderr}")
        print(f"Files after simulation: {os.listdir(work_dir)}")

        # 3. Convert trajectory to PDB
        traj_file = output_prefix.with_suffix(".trafl")
        if not traj_file.exists():
            raise FileNotFoundError(f"Trajectory file missing: {traj_file}")
        print(f"Found trajectory file: {traj_file}")

        # Convert with error capture
        print("Converting trajectory...")
        result_convert = subprocess.run([
            str(SIMRNA_TRAFL2PDBS),
            str(init_pdb),
            str(traj_file),
            "1"
        ], capture_output=True, text=True, check=True)
        
        print(f"Conversion stdout: {result_convert.stdout}")
        print(f"Conversion stderr: {result_convert.stderr}")

        # Handle output file
        output_frame = Path(f"{traj_file.stem}_1.pdb")
        if not output_frame.exists():
            # Try alternative naming patterns
            alt_output = next(Path().glob("*_1.pdb"), None)
            if not alt_output:
                raise FileNotFoundError("No converted PDB found")
            output_frame = alt_output
        
        print(f"Found converted PDB at: {output_frame}")
        shutil.copy(output_frame, output_pdb)
        print(f"Final PDB copied to: {output_pdb}")

    except subprocess.CalledProcessError as e:
        print(f"Subprocess error: {e}\nOutput: {e.output}\nStderr: {e.stderr}")
        generate_dummy_pdb(sequence, output_pdb)
    except Exception as e:
        print(f"General error: {str(e)}")
        generate_dummy_pdb(sequence, output_pdb)

    return {
        "secondary_structure": ss,
        "mfe": mfe,
        "model_path": str(output_pdb)
    }


def predict_structure_with_rosetta(sequence):
    """Wrapper for compatibility"""
    return predict_aptamer_structure(sequence)
