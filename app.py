# app.py

import streamlit as st
import pandas as pd
import numpy as np
import py3Dmol
import requests

st.set_page_config(layout="wide")
st.title("🔬 Redox Structure Viewer")

# === Upload Data ===
uploaded_file = st.file_uploader("Upload redox site file (.tsv)", type=["tsv"])

if uploaded_file:
    df = pd.read_csv(uploaded_file, sep='\t')
    required_cols = {'Protein', 'Site', 'Sample', '%Reduced'}
    
    if not required_cols.issubset(df.columns):
        st.error(f"Missing required columns. Found: {df.columns.tolist()}")
        st.stop()

    # === Preprocess ===
    df['%Oxidized'] = 100 - df['%Reduced']
    df['Condition'] = df['Sample'].str.extract(r'(Fresh|Store)', expand=False)
    df = df.dropna(subset=['Protein', 'Site', '%Oxidized', 'Condition'])

    # === Summarize ===
    summary = (
        df.groupby(['Protein', 'Site', 'Condition'])['%Oxidized']
        .agg(['mean', 'std'])
        .reset_index()
        .pivot(index=['Protein', 'Site'], columns='Condition')
        .reset_index()
    )

    # Flatten MultiIndex columns
    summary.columns = ['Protein', 'Site', 'Fresh_mean', 'Store_mean', 'Fresh_std', 'Store_std']
    summary['DeltaOx'] = summary['Fresh_mean'] - summary['Store_mean']

    # === Select Protein ===
    uniprot_list = sorted(summary['Protein'].unique())
    selected_protein = st.selectbox("Select UniProt ID", uniprot_list)

    protein_df = summary[summary['Protein'] == selected_protein]

    # === Plot summary ===
    st.subheader("📊 Oxidation Summary Per Site")
    st.dataframe(protein_df[['Site', 'Fresh_mean', 'Store_mean', 'DeltaOx']].round(2))

    # === AlphaFold fetch ===
    st.subheader("🧬 AlphaFold Structure")
    af_url = f"https://alphafold.ebi.ac.uk/files/AF-{selected_protein}-F1-model_v4.pdb"
    r = requests.get(af_url)

    if r.status_code != 200:
        st.error("AlphaFold structure not found.")
        st.stop()

    pdb_lines = r.text
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_lines, 'pdb')
    viewer.setStyle({'cartoon': {'color': 'white'}})

    # === Color by Oxidation ===
    st.markdown("### 🎨 Color By:")
    color_option = st.radio("Choose coloring metric:", ['Fresh_mean', 'Store_mean', 'DeltaOx'])

    cmap = plt.get_cmap("coolwarm")
    values = protein_df[color_option].fillna(0)
    vmin, vmax = values.min(), values.max()

    for _, row in protein_df.iterrows():
        site = int(row['Site'])
        val = row[color_option]
        if pd.notna(val):
            norm_val = (val - vmin) / (vmax - vmin + 1e-6)
            r, g, b, _ = cmap(norm_val)
            hex_color = '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))
            viewer.addStyle({'resi': site}, {'cartoon': {'color': hex_color}})

    viewer.zoomTo()
    viewer.show()
    st.components.v1.html(viewer.render(), height=600)
