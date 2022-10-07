"""
OpenFold Implementation in a combination of the colab notebook and the original
"""
import unittest.mock
import sys
import os

from urllib import request
from concurrent import futures
import json

# from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np

# import py3Dmol
from tqdm import tqdm
import torch

# from IPython import display
# from ipywidgets import GridspecLayout
# from ipywidgets import Output

from openfold import config
from openfold.data import feature_pipeline
from openfold.data import parsers
from openfold.data import data_pipeline
from openfold.data.tools import jackhmmer
from openfold.model import model
from openfold.np import protein
from openfold.np.relax import relax
from openfold.np.relax import utils
from openfold.utils.tensor_utils import tensor_tree_map


# A filthy hack to avoid slow Linear layer initialization
import openfold.model.primitives

from fastapi import FastAPI

app = FastAPI()

def __default_linear_init__(self, *args, **kwargs):
    return torch.nn.Linear.__init__(
        self, *args[:2], **{k: v for k, v in kwargs.items() if k == "bias"}
    )


openfold.model.primitives.Linear.__init__ = __default_linear_init__


RELAX_PREDICTION = True

UNEDITED_SEQUENCE = ""

with open("/openfold/input.txt", encoding="utf-8") as f:
    UNEDITED_SEQUENCE = f.readline()

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

sequence = format_input(UNEDITED_SEQUENCE)

@app.get("/sequence/{sequence}")
def read_sequence(sequence: str):
    return {"sequence": format_input(sequence)}

sys.path.insert(0, "/usr/local/lib/python3.7/site-packages/")
sys.path.append("/opt/conda/lib/python3.7/site-packages")

# Allows us to skip installing these packages
unnecessary_modules = [
    "dllogger",
    "pytorch_lightning",
    "pytorch_lightning.utilities",
    "pytorch_lightning.callbacks.early_stopping",
    "pytorch_lightning.utilities.seed",
]
for unnecessary_module in unnecessary_modules:
    sys.modules[unnecessary_module] = unittest.mock.MagicMock()

# @title Search against genetic databases

# @markdown Once this cell has been executed, you will see
# @markdown statistics about the multiple sequence alignment
# @markdown (MSA) that will be used by OpenFold. In particular,
# @markdown you’ll see how well each residue is covered by similar
# @markdown sequences in the MSA.

# --- Find the closest source ---
TEST_URL_PATTERN = (
    "https://storage.googleapis.com/alphafold-colab{:s}"
    + "/latest/uniref90_2021_03.fasta.1"
)
ex = futures.ThreadPoolExecutor(3)

def fetch(source):
    """--- Find the closest source ---"""
    request.urlretrieve(TEST_URL_PATTERN.format(source))
    return source


fs = [ex.submit(fetch, source) for source in ["", "-europe", "-asia"]]
SOURCE = None
for f in futures.as_completed(fs):
    SOURCE = f.result()
    ex.shutdown()
    break

# --- Search against genetic databases ---
with open("target.fasta", "wt", encoding="utf-8") as f:
    f.write(f">query\n{sequence}")

# Run the search against chunks of genetic databases (since the genetic
# databases don't fit in Colab ramdisk).

JACKHMMER_BINARY_PATH = "/usr/bin/jackhmmer"
dbs = []
TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt}"
    + " [elapsed: {elapsed} remaining: {remaining}]"
)

NUM_JACKHMMER_CHUNKS = {"uniref90": 59, "smallbfd": 17, "mgnify": 71}
TOTAL_JACKHMMER_CHUNKS = sum(NUM_JACKHMMER_CHUNKS.values())
with tqdm(total=TOTAL_JACKHMMER_CHUNKS, bar_format=TQDM_BAR_FORMAT) as pbar:

    def jackhmmer_chunk_callback(i):
        """Update the tqdm bar progress"""
        pbar.update(n=1)

    pbar.set_description("Searching uniref90")
    jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
        binary_path=JACKHMMER_BINARY_PATH,
        database_path=f"https://storage.googleapis.com/alphafold\
            -colab{SOURCE}/latest/uniref90_2021_03.fasta",
        get_tblout=True,
        num_streamed_chunks=NUM_JACKHMMER_CHUNKS["uniref90"],
        streaming_callback=jackhmmer_chunk_callback,
        z_value=135301051,
    )
    dbs.append(("uniref90", jackhmmer_uniref90_runner.query("target.fasta")))

    pbar.set_description("Searching smallbfd")
    jackhmmer_smallbfd_runner = jackhmmer.Jackhmmer(
        binary_path=JACKHMMER_BINARY_PATH,
        database_path=f"https://storage.googleapis.com/alphafold\
            -colab{SOURCE}/latest/bfd-first_non_consensus_sequences.fasta",
        get_tblout=True,
        num_streamed_chunks=NUM_JACKHMMER_CHUNKS["smallbfd"],
        streaming_callback=jackhmmer_chunk_callback,
        z_value=65984053,
    )
    dbs.append(("smallbfd", jackhmmer_smallbfd_runner.query("target.fasta")))

    pbar.set_description("Searching mgnify")
    jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
        binary_path=JACKHMMER_BINARY_PATH,
        database_path=f"https://storage.googleapis.com/alphafold\
            -colab{SOURCE}/latest/mgy_clusters_2019_05.fasta",
        get_tblout=True,
        num_streamed_chunks=NUM_JACKHMMER_CHUNKS["mgnify"],
        streaming_callback=jackhmmer_chunk_callback,
        z_value=304820129,
    )
    dbs.append(("mgnify", jackhmmer_mgnify_runner.query("target.fasta")))


# --- Extract the MSAs and visualize ---
# Extract the MSAs from the Stockholm files.
# NB: deduplication happens later in data_pipeline.make_msa_features.

MGNIFY_MAX_HITS = 501

msas = []
deletion_matrices = []
full_msa = []
for db_name, db_results in dbs:
    unsorted_results = []
    for i, result in enumerate(db_results):
        msa, deletion_matrix, target_names = parsers.parse_stockholm(
            result["sto"]
        )
        e_values_dict = parsers.parse_e_values_from_tblout(result["tbl"])
        e_values = [e_values_dict[t.split("/")[0]] for t in target_names]
        zipped_results = zip(msa, deletion_matrix, target_names, e_values)
        if i != 0:
            # Only take query from the first chunk
            zipped_results = [x for x in zipped_results if x[2] != "query"]
        unsorted_results.extend(zipped_results)
    sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[3])
    db_msas, db_deletion_matrices, _, _ = zip(*sorted_by_evalue)
    if db_msas:
        if db_name == "mgnify":
            db_msas = db_msas[:MGNIFY_MAX_HITS]
            db_deletion_matrices = db_deletion_matrices[:MGNIFY_MAX_HITS]
        full_msa.extend(db_msas)
        msas.append(db_msas)
        deletion_matrices.append(db_deletion_matrices)
        msa_size = len(set(db_msas))
        print(f"{msa_size} Sequences Found in {db_name}")

deduped_full_msa = list(dict.fromkeys(full_msa))
TOTAL_MSA_SIZE = len(deduped_full_msa)
print(f"\n{TOTAL_MSA_SIZE} Sequences Found in Total\n")

aa_map = {
    restype: i for i, restype in enumerate(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ-"
    )
}
msa_arr = np.array([[aa_map[aa] for aa in seq] for seq in deduped_full_msa])
num_alignments, num_res = msa_arr.shape

fig = plt.figure(figsize=(12, 3))
plt.title("Per-Residue Count of Non-Gap Amino Acids in the MSA")
plt.plot(np.sum(msa_arr != aa_map["-"], axis=0), color="black")
plt.ylabel("Non-Gap Count")
plt.yticks(range(0, num_alignments + 1, max(1, int(num_alignments / 3))))
plt.savefig("msa.png", dpi=300, bbox_inches="tight")

# @title Run OpenFold and download prediction

# @markdown Once this cell has been executed, a zip-archive with
# @markdown the obtained prediction will be automatically downloaded
# @markdown to your computer.

# Color bands for visualizing plddt
PLDDT_BANDS = [
    (0, 50, "#FF7D45"),
    (50, 70, "#FFDB13"),
    (70, 90, "#65CBF3"),
    (90, 100, "#0053D6"),
]

# --- Run the model ---
model_names = [
    "finetuning_3.pt",
    "finetuning_4.pt",
    "finetuning_5.pt",
    "finetuning_ptm_2.pt",
    "finetuning_no_templ_ptm_1.pt",
]


def _placeholder_template_feats(num_templates_, num_res_):
    return {
        "template_aatype": np.zeros(
            (num_templates_, num_res_, 22), dtype=np.int64
        ),
        "template_all_atom_positions": np.zeros(
            (num_templates_, num_res_, 37, 3), dtype=np.float32
        ),
        "template_all_atom_mask": np.zeros(
            (num_templates_, num_res_, 37), dtype=np.float32
        ),
        "template_domain_names": np.zeros((num_templates_,), dtype=np.float32),
        "template_sum_probs": np.zeros((num_templates_, 1), dtype=np.float32),
    }


OUTPUT_DIR = "prediction"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plddts = {}
pae_outputs = {}
unrelaxed_proteins = {}
OPENFOLD_PARAMS_DIR = "/content/openfold/openfold/resources/openfold_params"

with tqdm(total=len(model_names) + 1, bar_format=TQDM_BAR_FORMAT) as pbar:
    for i, model_name in list(enumerate(model_names)):
        pbar.set_description(f"Running {model_name}")
        NUM_TEMPLATES = 1  # dummy number --- is ignored
        num_res = len(sequence)

        feature_dict = {}
        feature_dict.update(
            data_pipeline.make_sequence_features(
                sequence,
                "test",
                num_res
            )
        )
        feature_dict.update(
            data_pipeline.make_msa_features(
                msas,
                deletion_matrices=deletion_matrices
            )
        )
        feature_dict.update(
            _placeholder_template_feats(NUM_TEMPLATES, num_res)
        )

        if "_no_templ_" in model_name:
            CONFIG_PRESET = "model_3"
        else:
            CONFIG_PRESET = "model_1"
        if "_ptm_" in model_name:
            CONFIG_PRESET += "_ptm"

        cfg = config.model_config(CONFIG_PRESET)
        openfold_model = model.AlphaFold(cfg)
        openfold_model = openfold_model.eval()
        params_name = os.path.join(OPENFOLD_PARAMS_DIR, model_name)
        d = torch.load(params_name)
        openfold_model.load_state_dict(d)

        openfold_model = openfold_model.cuda()

        pipeline = feature_pipeline.FeaturePipeline(cfg.data)
        processed_feature_dict = pipeline.process_features(
            feature_dict,
            mode="predict"
        )

        processed_feature_dict = tensor_tree_map(
            lambda t: t.cuda(), processed_feature_dict
        )

        with torch.no_grad():
            prediction_result = openfold_model(processed_feature_dict)

        # Move the batch and output to np for further processing
        processed_feature_dict = tensor_tree_map(
            lambda t: np.array(t[..., -1].cpu()), processed_feature_dict
        )
        prediction_result = tensor_tree_map(
            lambda t: np.array(t.cpu()), prediction_result
        )

        mean_plddt = prediction_result["plddt"].mean()

        if "predicted_aligned_error" in prediction_result:
            pae_outputs[model_name] = (
                prediction_result["predicted_aligned_error"],
                prediction_result["max_predicted_aligned_error"],
            )
        else:
            # Get the pLDDT confidence metrics.
            # Do not put pTM models here as they
            # should never get selected.
            plddts[model_name] = prediction_result["plddt"]

        # Set the b-factors to the per-residue plddt.
        final_atom_mask = prediction_result["final_atom_mask"]
        b_factors = prediction_result["plddt"][:, None] * final_atom_mask
        unrelaxed_protein = protein.from_prediction(
            processed_feature_dict, prediction_result, b_factors=b_factors
        )
        unrelaxed_proteins[model_name] = unrelaxed_protein

        # Delete unused outputs to save memory.
        del openfold_model
        del processed_feature_dict
        del prediction_result
        pbar.update(n=1)

    # Find the best model according to the mean pLDDT.
    best_model_name = max(plddts.keys(), key=lambda x: plddts[x].mean())
    BEST_PDB = protein.to_pdb(unrelaxed_proteins[best_model_name])

    # --- AMBER relax the best model ---
    if RELAX_PREDICTION:
        pbar.set_description(f"AMBER relaxation")
        amber_relaxer = relax.AmberRelaxation(
            max_iterations=0,
            tolerance=2.39,
            stiffness=10.0,
            exclude_residues=[],
            max_outer_iterations=20,
            use_gpu=False,
        )
        relaxed_pdb, _, _ = amber_relaxer.process(
            prot=unrelaxed_proteins[best_model_name]
        )

        # Write out the prediction
        pred_output_path = os.path.join(OUTPUT_DIR, "selected_prediction.pdb")
        with open(pred_output_path, "w", encoding="utf-8") as f:
            f.write(relaxed_pdb)

        BEST_PDB = relaxed_pdb

    pbar.update(n=1)  # Finished AMBER relax.

# Construct multiclass b-factors to indicate confidence bands
# 0=very low, 1=low, 2=confident, 3=very high
banded_b_factors = []
for plddt in plddts[best_model_name]:
    for idx, (min_val, max_val, _) in enumerate(PLDDT_BANDS):
        if min_val <= plddt <= max_val:
            banded_b_factors.append(idx)
            break
banded_b_factors = np.array(banded_b_factors)[:, None] * final_atom_mask
to_visualize_pdb = utils.overwrite_b_factors(BEST_PDB, banded_b_factors)

# --- Visualise the prediction & confidence ---
SHOW_SIDECHAINS = True
# def plot_plddt_legend():
#     """Plots the legend for pLDDT."""
#     thresh = [
#         'Very low (pLDDT < 50)',
#         'Low (70 > pLDDT > 50)',
#         'Confident (90 > pLDDT > 70)',
#         'Very high (pLDDT > 90)']

#     colors = [x[2] for x in PLDDT_BANDS]

#     plt.figure(figsize=(2, 2))
#     for c in colors:
#         plt.bar(0, 0, color=c)
#     plt.legend(thresh, frameon=False, loc='center', fontsize=20)
#     plt.xticks([])
#     plt.yticks([])
#     ax = plt.gca()
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.spines['left'].set_visible(False)
#     ax.spines['bottom'].set_visible(False)
#     plt.title('Model Confidence', fontsize=20, pad=20)
#     return plt


# Visualise the prediction here
# Doesn't really work with the terminal in docker

# Color the structure by per-residue pLDDT
# color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}
# view = py3Dmol.view(width=800, height=600)
# view.addModelsAsFrames(to_visualize_pdb)
# style = {'cartoon': {
#     'colorscheme': {
#         'prop': 'b',
#         'map': color_map}
#         }}
# if SHOW_SIDECHAINS:
#     style['stick'] = {}
# view.setStyle({'model': -1}, style)
# view.zoomTo()

# grid = GridspecLayout(1, 2)
# out = Output()
# with out:
#     view.savefig('prediction.png')
#     view.show()
# grid[0, 0] = out

# out = Output()
# with out:
#     plot_plddt_legend().show()
# grid[0, 1] = out

# display.display(grid)

# Display pLDDT and predicted aligned error (if output by the model).
if pae_outputs:
    NUM_PLOTS = 2
else:
    NUM_PLOTS = 1

plt.figure(figsize=[8 * NUM_PLOTS, 6])
plt.subplot(1, NUM_PLOTS, 1)
plt.plot(plddts[best_model_name])
plt.title("Predicted LDDT")
plt.xlabel("Residue")
plt.ylabel("pLDDT")

if NUM_PLOTS == 2:
    plt.subplot(1, 2, 2)
    pae, max_pae = list(pae_outputs.values())[0]
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title("Predicted Aligned Error")
    plt.xlabel("Scored residue")
    plt.ylabel("Aligned residue")

plt.imsave("prediction.png", pae, vmin=0.0, vmax=max_pae, cmap="Greens_r")

# Save pLDDT and predicted aligned error (if it exists)
pae_output_path = os.path.join(OUTPUT_DIR, "predicted_aligned_error.json")
if pae_outputs:
    # Save predicted aligned error in the same format as the AF EMBL DB
    rounded_errors = np.round(pae.astype(np.float64), decimals=1)
    indices = np.indices((len(rounded_errors), len(rounded_errors))) + 1
    indices_1 = indices[0].flatten().tolist()
    indices_2 = indices[1].flatten().tolist()
    pae_data = json.dumps(
        [
            {
                "residue1": indices_1,
                "residue2": indices_2,
                "distance": rounded_errors.flatten().tolist(),
                "max_predicted_aligned_error": max_pae.item(),
            }
        ],
        indent=None,
        separators=(",", ":"),
    )
    with open(pae_output_path, "w", encoding="utf-8") as f:
        f.write(pae_data)
