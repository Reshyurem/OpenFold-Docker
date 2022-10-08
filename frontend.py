import streamlit as st
import py3Dmol
from stmol import showmol
st.sidebar.title('Show Proteins')
prot_str='1A2C,1BML,1D5M,1D5X,1D5Z,1D6E,1DEE,1E9F,1FC2,1FCC,1G4U,1GZS,1HE1,1HEZ,1HQR,1HXY,1IBX,1JBU,1JWM,1JWS'
prot_list=prot_str.split(',')
protein = st.sidebar.text_input('Enter a protein name', '1A2C')
spin = st.sidebar.checkbox('Spin', value = False)
xyzview = py3Dmol.view(query='pdb:'+protein)
xyzview.setStyle({'cartoon':{'color':'spectrum'}})
xyzview.setBackgroundColor('#FFFFFF')
if spin:
    xyzview.spin(True)
else:
    xyzview.spin(False)
xyzview.zoomTo()
showmol(xyzview,height=500,width=800)