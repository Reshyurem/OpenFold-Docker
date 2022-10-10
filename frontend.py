import streamlit as st
import py3Dmol
from stmol import showmol
st.sidebar.title('Openfold')
ml = st.sidebar.form("Protein Model Generation")
protein = ml.text_input('Enter a protein sequence')
generate = ml.form_submit_button("Generate Protein PDB")

if generate:
    print("Hello")

def render_mol(pdb):
    pdbview = py3Dmol.view(width=400,height=400)
    pdbview.addModelsAsFrames(pdb)
    pdbview.setStyle({'cartoon':{'color':'spectrum'}})
    # pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    if spin:
        pdbview.spin(True)
    else:
        pdbview.spin(False)
    showmol(pdbview, height = 500,width=800)

spin = st.sidebar.checkbox('Spin', value = False)
uploaded_file = st.sidebar.file_uploader("Choose the file to model")
if uploaded_file is not None:
    # pdb = uploaded_file.getvalue().decode("utf-8")
    with open("cocaine.pdb") as ifile:
        pdb = "".join([x for x in ifile])

    render_mol(pdb)