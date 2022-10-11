import streamlit as st
import py3Dmol
from stmol import showmol
st.sidebar.title('Openfold')
ml = st.sidebar.form("Protein Model Generation")
protein = ml.text_input('Enter a protein sequence')
generate = ml.form_submit_button("Generate Protein PDB")

if generate:
    print("Hello")

# Color bands for visualizing plddt
PLDDT_BANDS = [
    (0, 50, "#FF7D45"),
    (50, 70, "#FFDB13"),
    (70, 90, "#65CBF3"),
    (90, 100, "#0053D6"),
]
color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}

def render_mol(pdb):
    pdbview = py3Dmol.view(width=400,height=400)
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