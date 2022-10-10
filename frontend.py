import streamlit as st
import py3Dmol
from stmol import showmol
st.sidebar.title('Openfold')
protein = st.sidebar.text_input('Enter a protein name', '1A2C')
spin = st.sidebar.checkbox('Spin', value = False)

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

uploaded_files = st.sidebar.file_uploader("Choose xyz files")
for uploaded_file in uploaded_files:
    xyz = uploaded_file.getvalue().decode("utf-8")
    render_mol(xyz)