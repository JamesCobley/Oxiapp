# app.py

import streamlit as st
import pandas as pd
import numpy as np
import py3Dmol
import requests

# === SETTINGS ===
RAW_DATA_URL = "https://raw.githubusercontent.com/JamesCobley/Oxiapp/main/site_redox_quant.tsv"

st.set_page_config(layout="wide")
st.title("🔬 Oxidation Viewer via AlphaFold")

@st.cache_data
def load_data(url):
    df = pd.read_csv(url, sep='\t')
    df['Site'] = df['Site'].astype(str).str.extract(r'(\d+)').astype(float).astype(int)
    return df

@st.cache_data
def fetch_alphafold_pdb(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.text

def color_by_oxidation(pdb, site_values, color_map='redox'):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb, 'pdb')
    view.setStyle({'cartoon': {'color': 'white'}})

    for site, value in site_values.items():
        if value is None: continue
        color = get_color(value)
        view.addStyle(
            {'resi': int(site), 'atom': 'SG'},  # Cys SG atom
            {'stick': {'color': color, 'radius': 0.3}}
        )
    view.zoomTo()
    return view

def get_color(percent_ox):
    if percent_ox < 20: return 'blue'
    elif percent_ox < 40: return 'green'
    elif percent_ox < 60: return 'yellow'
    elif percent_ox < 80: return 'orange'
    else: return 'red'

# === MAIN ===

df = load_data(RAW_DATA_URL)

uniprots = sorted(df['Protein'].dropna().unique())
choice = st.selectbox("Select UniProt ID", uniprots)

view_mode = st.radio("View", ['Fresh', 'Store', 'Delta'], horizontal=True)

# Subset data
subset = df[df['Protein'] == choice]

# Pick relevant samples
fresh_cols = [col for col in df.columns if 'Fresh' in col and col.endswith('%Reduced')]
store_cols = [col for col in df.columns if 'Store' in col and col.endswith('%Reduced')]

# Compute %Oxidized
subset['Fresh_Ox'] = 100 - subset[fresh_cols].mean(axis=1)
subset['Store_Ox'] = 100 - subset[store_cols].mean(axis=1)
subset['Delta'] = subset['Fresh_Ox'] - subset['Store_Ox']

# Choose data for coloring
if view_mode == 'Fresh':
    site_color_data = dict(zip(subset['Site'], subset['Fresh_Ox']))
elif view_mode == 'Store':
    site_color_data = dict(zip(subset['Site'], subset['Store_Ox']))
else:
    site_color_data = dict(zip(subset['Site'], subset['Delta']))

# Fetch structure
with st.spinner("Fetching AlphaFold structure..."):
    pdb = fetch_alphafold_pdb(choice)

if pdb:
    view = color_by_oxidation(pdb, site_color_data)
    view.show()
    view.png()
    st.pydeck_chart(view)
else:
    st.error("❌ AlphaFold structure not found for this UniProt ID.")

# Show data table
st.subheader("📋 Site-level Data")
st.dataframe(subset[['Site', 'Fresh_Ox', 'Store_Ox', 'Delta']])

