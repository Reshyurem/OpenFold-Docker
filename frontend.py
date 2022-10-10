import streamlit as st
import py3Dmol
from stmol import showmol
st.sidebar.title('Openfold')
ml = st.sidebar.form("Protein Model Generation")
protein = ml.text_input('Enter a protein sequence')
generate = ml.form_submit_button("Generate Protein PDB")

if generate:
    print("Hello")

def render_mol(xyz):
    xyzview = py3Dmol.view(width=400,height=400)
    xyzview.addModel(xyz,'xyz')
    xyzview.setStyle({'color':'spectrum'})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    if spin:
        xyzview.spin(True)
    else:
        xyzview.spin(False)
    showmol(xyzview, height = 500,width=800)

model = st.sidebar.container
spin = st.sidebar.checkbox('Spin', value = False)
uploaded_file = st.sidebar.file_uploader("Choose the file to model")
if uploaded_file is not None:
    xyz = uploaded_file.getvalue().decode("utf-8")
    render_mol(xyz)