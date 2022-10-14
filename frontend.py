import streamlit as st
import requests
import py3Dmol
import time
from stmol import showmol

def format_input(unformatted_input):
    """Function to format the input sequence"""
    formatted_input = unformatted_input.translate(str.maketrans("", "", " \n\t")).upper()
    aatypes = set("ACDEFGHIKLMNPQRSTVWY")  # 20 standard aatypes

    if not set(formatted_input).issubset(aatypes):
        raise Exception(
            f"Input sequence contains non-amino acid letters: \
            {set(formatted_input) - aatypes}. \
            OpenFold only supports 20 standard amino acids as inputs."
        )

    return formatted_input

st.sidebar.title('Openfold')
ml = st.sidebar.form("Protein Model Generation")
protein = format_input(ml.text_input('Enter a protein sequence'))
generate = ml.form_submit_button("Generate Protein PDB")

if generate and (protein is not None):
    r = requests.get(url = "http://localhost:8000/sequence/" + protein)
    st.sidebar.download_button('Click to download the file', r.content, file_name='protein.pdb')

# Color bands for visualizing plddt
PLDDT_BANDS = [
    (0, 50, "#FF7D45"),
    (50, 70, "#FFDB13"),
    (70, 90, "#65CBF3"),
    (90, 100, "#0053D6"),
]
color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}

def render_mol(pdb):
    pdbview = py3Dmol.view(width=600,height=600)
    pdbview.addModelsAsFrames(pdb)
    style = {'cartoon': {
        'colorscheme': {
            'prop': 'b',
            'map': color_map}
            }}
    style['stick'] = {}
    pdbview.setStyle(style)
    pdbview.zoomTo()
    if spin:
        pdbview.spin(True)
    else:
        pdbview.spin(False)
    showmol(pdbview, height = 500,width=800)

spin = st.sidebar.checkbox('Spin', value = False)
uploaded_file = st.sidebar.file_uploader("Choose the file to model")
if uploaded_file is not None:
    pdb = uploaded_file.getvalue().decode("utf-8")
    render_mol(pdb)