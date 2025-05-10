import streamlit as st
import pandas as pd
import subprocess
from pathlib import Path
from Bio.SeqUtils import molecular_weight
from backend.generate import generate_initial_candidates
from backend.structure import predict_secondary_structure
from backend.docking import run_docking_analysis
from backend.optimize import optimize_candidates

# Configuration
TEMP_DIR = Path("temp")
TEMP_DIR.mkdir(exist_ok=True)

st.set_page_config(page_title="Aptamer AI Designer", layout="wide")

def main():
    st.title("Open Source AI Aptamer Designer")
    st.markdown("### ecoSPECS x GC.OS Collaboration Platform")
    
    # Sidebar for inputs
    with st.sidebar:
        st.header("Configuration")
        uploaded_file = st.file_uploader("Upload Target Protein (PDB)", type=["pdb"])
        seq_length = st.slider("Aptamer Length", 15, 40, 20)
        num_candidates = st.number_input("Number of Candidates", 5, 100, 20)
        run_docking = st.checkbox("Enable Docking Analysis", True)
    
    # Main interface
    if uploaded_file and st.button("Start Design Pipeline"):
        target_path = TEMP_DIR / uploaded_file.name
        target_path.write_bytes(uploaded_file.getvalue())
        
        with st.status("Running Aptamer Design Pipeline...", expanded=True) as status:
            # Stage 1: Candidate Generation
            st.write("1. Generating initial candidates...")
            candidates = generate_initial_candidates(num_candidates, seq_length)
            
            # Stage 2: Structure Prediction
            st.write("2. Predicting secondary structures...")
            structures = []
            for seq in candidates:
                try:
                    ss, mfe = predict_secondary_structure(seq)
                    structures.append({
                        "sequence": seq,
                        "structure": ss,
                        "mfe": mfe,
                        "weight": molecular_weight(seq, seq_type="DNA")
                    })
                except Exception as e:
                    st.error(f"Error processing {seq}: {str(e)}")
            df = pd.DataFrame(structures)
            
            # Stage 3: Docking Analysis
            if run_docking:
                st.write("3. Running docking analysis...")
                try:
                    docking_results = run_docking_analysis(df, target_path)
                    df = pd.merge(df, docking_results, on="sequence")
                except Exception as e:
                    st.error(f"Docking failed: {str(e)}")
            
            # Stage 4: Optimization
            st.write("4. Optimizing candidates...")
            optimized = optimize_candidates(df)
            status.update(label="Pipeline complete!", state="complete", expanded=False)
        
        # Display results
        st.subheader("Top Candidates")
        cols = st.columns(3)
        cols[0].dataframe(optimized.sort_values("score", ascending=False).head(5))
        
        with cols[1]:
            st.plotly_chart(plot_energy_distribution(optimized))
        
        with cols[2]:
            st.download_button(
                label="Download Results",
                data=optimized.to_csv(index=False),
                file_name="aptamer_candidates.csv"
            )

def plot_energy_distribution(df):
    import plotly.express as px
    fig = px.histogram(df, x="mfe", nbins=20, 
                      title="Free Energy Distribution")
    fig.update_layout(height=300)
    return fig

if __name__ == "__main__":
    main()
