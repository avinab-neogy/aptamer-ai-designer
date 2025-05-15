import streamlit as st
import sys
import os
import traceback
import base64
from pathlib import Path
import pandas as pd
import subprocess
from Bio.SeqUtils import molecular_weight
import hashlib
import plotly.express as px

# Set background image
def set_background(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode()
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/jpg;base64,{encoded_string}");
            background-size: cover;
            background-attachment: fixed;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

st.set_page_config(
    page_title="Aptamer AI Designer", 
    layout="wide",
    menu_items={
        'Get Help': 'https://github.com/your-repo',
        'Report a bug': "https://github.com/your-repo/issues",
        'About': "# Aptamer AI Designer\nAn open-source aptamer design platform"
    }
)

# Apply core UI styling
st.markdown("""
    <style>
    .block-container {
        background: rgba(255, 255, 255, 0.4);
        border-radius: 10px;
        padding: 2rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        backdrop-filter: blur(8px);
        margin-top: 2rem;
    }
    .stButton>button {
        background: linear-gradient(45deg, #6366f1 0%, #60a5fa 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.75rem 1.5rem;
        font-weight: 600;
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    .stSidebar {
        background: rgba(255, 255, 255, 0.85) !important;
        backdrop-filter: blur(6px);
        border-right: 1px solid rgba(0, 0, 0, 0.1);
    }
    [data-testid=stSidebar] * {
        color: #222 !important;
    }
    .stSlider > div, .stNumberInput > div, .stCheckbox > label, .stFileUploader > div {
        color: #222 !important;
    }
    .stMarkdown h2, .stMarkdown h3, .stMarkdown h4, .stMarkdown h5, .stMarkdown h6 {
        color: #222 !important;
    }
    .stSlider > div [role=slider] {
        background: #222 !important;
        border: 2px solid #333 !important;
    }
    h1, h3 {
        color: #222 !important;
    }
    </style>
""", unsafe_allow_html=True)

# --- Enhanced file uploader styling for light box ---
st.markdown("""
<style>
[data-testid="stFileUploadDropzone"] {
    background: #fff !important;
    border: 2px solid #60a5fa !important;
    border-radius: 8px !important;
    box-shadow: 0 2px 12px rgba(0,0,0,0.06);
    padding-top: 10px !important;
    padding-bottom: 10px !important;
}
[data-testid="stFileUploadDropzone"] * {
    color: #222 !important;
    font-weight: 600 !important;
    text-shadow: none !important;
}
[data-testid="stFileUploadDropzone"] button {
    color: #fff !important;
    background: #6366f1 !important;
    border-radius: 6px !important;
    font-weight: 700 !important;
    border: none !important;
    box-shadow: 0 1px 6px rgba(99,102,241,0.10);
}
</style>
""", unsafe_allow_html=True)

set_background("assets/background.jpg")

# Absolute path to project root
PROJECT_ROOT = Path("/home/avinab/Documents/aptamer_ai")
sys.path.append(str(PROJECT_ROOT))

from backend.generate import generate_initial_candidates
from backend.structure import predict_secondary_structure, predict_structure_with_rosetta
from backend.docking import run_docking_analysis
from backend.optimize import optimize_candidates

# Configuration
TEMP_DIR = Path("temp")
TEMP_DIR.mkdir(exist_ok=True, parents=True)
os.chmod(TEMP_DIR, 0o755)

def render_svg(svg_path, width=140):
    with open(svg_path, "r") as f:
        svg = f.read()
    b64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
    html = f'<img src="data:image/svg+xml;base64,{b64}" width="{width}"/>'
    st.markdown(html, unsafe_allow_html=True)

def main():
    st.markdown('<h1 style="color:#222;font-weight:900;">🧬 Aptamer AI Designer</h1>', unsafe_allow_html=True)
    st.markdown('<h3 style="color:#333;font-weight:700;">Accelerated Aptamer Discovery Platform</h3>', unsafe_allow_html=True)
    
    # Sidebar for inputs
    with st.sidebar:
        logo_path = Path("assets/logo.svg")
        try:
            render_svg(logo_path, width=140)
        except Exception as e:
            st.warning(f"ecoSPECS logo could not be loaded: {e}")
        st.header("⚙️ Configuration")
        uploaded_file = st.file_uploader("Upload Target Protein (PDB)", type=["pdb"],
                                        help="Upload your target protein structure in PDB format")
        seq_length = st.slider("Aptamer Length", 15, 40, 20,
                                help="Select desired length for generated aptamers")
        num_candidates = st.number_input("Number of Candidates", 5, 100, 20,
                                        help="Number of initial candidates to generate")
        run_docking = st.checkbox("Enable Docking Analysis", True,
                                help="Perform molecular docking analysis")
        enable_3d = st.checkbox("Enable 3D Structure Prediction", False,
                                help="Generate 3D structural models")
    
    # Main interface
    if uploaded_file and st.button("🚀 Start Design Pipeline"):
        target_path = TEMP_DIR / uploaded_file.name
        target_path.write_bytes(uploaded_file.getvalue())
        
        with st.status("🔍 Running Aptamer Design Pipeline...", expanded=True) as status:
            # Stage 1: Candidate Generation
            st.write("1. 🧪 Generating initial candidates...")
            candidates = generate_initial_candidates(num_candidates, seq_length)
            
            # Stage 2: Structure Prediction
            st.write("2. 🧬 Predicting structures...")
            structures = []
            for seq in candidates:
                try:
                    if enable_3d:
                        result = predict_structure_with_rosetta(seq)
                        structures.append({
                            "sequence": seq,
                            "structure": result["secondary_structure"],
                            "mfe": result["mfe"],
                            "weight": molecular_weight(seq, seq_type="DNA"),
                            "model_path": result["model_path"]
                        })
                    else:
                        ss, mfe = predict_secondary_structure(seq)
                        structures.append({
                            "sequence": seq,
                            "structure": ss,
                            "mfe": mfe,
                            "weight": molecular_weight(seq, seq_type="DNA")
                        })
                except Exception as e:
                    st.error(f"Error processing {seq}: {str(e)}")
                    st.code(f"Error details:\n{traceback.format_exc()}", language='text')
                    structures.append({
                        "sequence": seq,
                        "structure": "",
                        "mfe": 0.0,
                        "weight": 0.0,
                        "model_path": None
                    })
            
            # Create DataFrame with fallback columns
            df = pd.DataFrame(structures) if structures else pd.DataFrame(columns=["sequence", "structure", "mfe", "weight", "model_path"])
            
            # Early exit if no valid candidates
            if df.empty:
                st.error("❌ No valid candidates generated. Please check your inputs.")
                return
            
            # Stage 3: Docking Analysis
            if run_docking:
                st.write("3. 🤖 Running docking analysis...")
                try:
                    docking_results = run_docking_analysis(df, target_path)
                    df = pd.merge(df, docking_results, on="sequence")
                except Exception as e:
                    st.error(f"Docking failed: {str(e)}")
            
            # Stage 4: Optimization
            st.write("4. 🚄 Optimizing candidates...")
            try:
                optimized = optimize_candidates(df)
            except KeyError as e:
                st.error(f"Optimization failed due to missing data: {str(e)}")
                optimized = df
            
            status.update(label="✅ Pipeline complete!", state="complete", expanded=False)
        
        # Display results
        st.subheader("🏆 Top Candidates")
        
        if not optimized.empty:
            cols = st.columns([2,3,2])
            with cols[0]:
                st.markdown("**Top Performers**")
                st.dataframe(optimized.sort_values("score", ascending=False).head(5),
                            use_container_width=True,
                            height=400)
            
            with cols[1]:
                st.markdown("**Energy Distribution**")
                st.plotly_chart(plot_energy_distribution(optimized),
                              use_container_width=True)
            
            with cols[2]:
                st.markdown("**Export Results**")
                st.download_button(
                    label="📥 Download CSV",
                    data=optimized.to_csv(index=False),
                    file_name="aptamer_candidates.csv",
                    mime="text/csv"
                )
                if enable_3d and 'model_path' in optimized.columns:
                    st.markdown("---")
                    st.markdown("**3D Visualization**")
                    selected_seq = st.selectbox(
                        "Select sequence", 
                        optimized.sort_values("score", ascending=False)["sequence"].tolist()
                    )
                    model_path = optimized[optimized["sequence"] == selected_seq]["model_path"].values[0]
                    
                    if model_path and Path(model_path).exists():
                        try:
                            import py3Dmol
                            view = py3Dmol.view(width=300, height=300)
                            with open(model_path, 'r') as f:
                                pdb_data = f.read()
                            view.addModel(pdb_data, 'pdb')
                            view.setStyle({'cartoon': {'color': 'spectrum'}})
                            view.zoomTo()
                            view.setBackgroundColor('0xeeeeee')
                            html_data = view._make_html()
                            st.components.v1.html(html_data, height=600)
                        except Exception as e:
                            st.error(f"Visualization error: {str(e)}")
                    else:
                        st.warning("3D model file not found")
        else:
            st.warning("⚠️ No optimized candidates to display")
        
        # Simulation diagnostics
        if enable_3d:
            with st.expander("🔬 Simulation Diagnostics"):
                for seq in candidates[:3]:
                    session_id = hashlib.md5(seq.encode()).hexdigest()[:8]
                    work_dir = TEMP_DIR/f"simrna_{session_id}"
                    if work_dir.exists():
                        st.markdown(f"**Files for `{seq[:10]}...`:**")
                        st.code("\n".join(os.listdir(work_dir)), language='text')
                    else:
                        st.warning(f"No directory found for {seq[:10]}...")

def plot_energy_distribution(df):
    fig = px.histogram(df, x="mfe", nbins=20, 
                      title="Free Energy Distribution",
                      color_discrete_sequence=['#6366f1'])
    fig.update_layout(
        height=300,
        plot_bgcolor='rgba(255,255,255,0.9)',
        paper_bgcolor='rgba(255,255,255,0.9)',
    )
    return fig

if __name__ == "__main__":
    main()
