#!/usr/bin/env python

import os
import io
import base64
import subprocess
from urllib.parse import quote
import tempfile
import statistics
import sys
import json


import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_table
# import dash_bio as dashbio
import plotly
import plotly.express as px
import plotly.graph_objects as go

import pandas as pd
from math import log

import Bio.SeqIO

from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser


from variant_data import (
    apply_mutations,
    graph_to_df,
    df_to_graph,
    SequenceInconsistency,
)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#              local app for exploration of protein engineering data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# start dashboard with DATA_PATH="example_data/gssm" python dashboard_local.py

###########################
# Data Manipulation / Model
###########################

# access environment variable, set default
DATA_PATH = os.getenv('DATA_PATH', os.path.join(os.path.dirname(os.path.realpath(__file__)), "example_data/directed"))
# print("path to data is", DATA_PATH)

# check if data path exists
if not os.path.isdir(DATA_PATH):
    sys.stderr.write("Error: data path does not exist\n")
    sys.exit(1)

CSV_FILE = os.path.join(DATA_PATH, "data.csv")
FASTA_FILE = os.path.join(DATA_PATH, "protein.fasta")
PDB_FILE = os.path.join(DATA_PATH, 'protein.pdb')
SETTINGS_FILE = os.path.join(DATA_PATH, 'settings.json')
HHBLITS_FASTA_FILE = os.path.join(DATA_PATH, "msa.fasta")
# HHBLITS_A3M_FILE = os.path.join(DATA_PATH, "4WOR_A_HHblits.a3m")

# check if all files are readible
for file_type in [CSV_FILE, FASTA_FILE, PDB_FILE, SETTINGS_FILE, HHBLITS_FASTA_FILE]:
    try:
        o = open(file_type,"r")
    except:
        sys.stderr.write("Error: Cannot open file: {}\n".format(file_type))
        sys.exit(1)

# load settings
SETTINGS = json.load(open(SETTINGS_FILE, "r"))
ATTRIBUTE_COL = SETTINGS["column_definitions"]["attribute"]
VARIANT_COL = SETTINGS["column_definitions"]["variant"]
PARENT_COL = SETTINGS["column_definitions"]["parent"]
MUTATION_COL = SETTINGS["column_definitions"]["mutation"]
# print(settings)

# OR keep track of headers but copy relevant columns, to new df with changed headers
# drop columns that are not specified? Would eventually 'clean' table

# read in sequence
for seqRec in Bio.SeqIO.parse(FASTA_FILE, "fasta"):
    wtSequence = seqRec.seq

# read in dataframe
df = pd.read_csv(CSV_FILE)
df = df.where(pd.notnull(df), None)
# drop columns that are not in settings
df = df[[VARIANT_COL, PARENT_COL, MUTATION_COL, ATTRIBUTE_COL]]

# rename columns the way graph functions expect it, attribute column won't be affected
df = df.rename(columns={
    VARIANT_COL: "variant",
    PARENT_COL: "parent",
    MUTATION_COL: "mutation"
})

try:
    print("df_to_graph")
    # returns set of variant objects
    variantGraph = df_to_graph(df, wtSequence)
    print("graph_to_df")
    # returns an extended version of the variant df
    df_extended = graph_to_df(df, variantGraph)
except SequenceInconsistency as e:
    print(e)


# rename columns back
df = df.rename(columns={
    "variant": VARIANT_COL,
    "parent": PARENT_COL,
    "mutation": MUTATION_COL
})


################################################
#               conservation plot
################################################

# subprocess.call(["./reformat.pl", 'a3m', 'fas', HHBLITS_A3M_FILE, HHBLITS_FASTA_FILE, '-r'])

# calculate information content for every position in the alignment 

# turn alignment into df
seq_df = pd.DataFrame([list(seqRec.seq._data) for seqRec in Bio.SeqIO.parse(HHBLITS_FASTA_FILE, "fasta")])
seq_length = len(seq_df.iloc[0].values)
if seq_length != len(wtSequence):
    raise Exception('homologous sequence and wt sequence have different length')

num_sequences = len(seq_df)
    
aa_alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

info_content = []
gap_frequency = []

# calculate information content according to https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation
samplesize_correction = (1 / log(2)) * ((20 - 1) / (2 * num_sequences))
# seq_len_normalization = 0
for i in range(seq_length):
    Hi = 0
    gap_frequency_site = seq_df.iloc[:,i].tolist().count('-') / num_sequences
    gap_frequency.append(gap_frequency_site * 100)
    for aa in aa_alphabet:
        freq = seq_df.iloc[:,i].tolist().count(aa) / (num_sequences)
        # print(freq)
        if freq != 0:
            Hi += freq * log(freq, 2)
    # account for gaps by multiplying with gep frequency
    info_content.append((log(20, 2) - (- Hi + samplesize_correction)) * (1 - gap_frequency_site))

# make df for plot
conservation_df = pd.DataFrame({'site': list(range(1, seq_length + 1)), 'conservation': info_content, 'gaps': gap_frequency})


def calculate_homology_logo():
    with tempfile.NamedTemporaryFile() as temp_file:
        homology_seqs = ','.join([''.join(row.tolist()) for i, row in seq_df.iterrows()])
        subprocess.call(["./ggseqlogo.R", temp_file.name, homology_seqs])
        data = base64.b64encode(open(temp_file.name, "rb").read())
        data_url_ggseqlogo_homology = "data:image/png;base64,{}".format(quote(data))
        return data_url_ggseqlogo_homology



####################################################
#           relative solvent accessibility      
####################################################
# http://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-module.html 

def check_for_shift(sequence, aa):
    '''accepts wtSequence, list of aa, list of positions'''
    # try every possible shift
    all_match = True
    for shift in range(len(sequence) - len(aa)):
        all_match = True
        # print('shift', shift)
        for j in range(len(aa) - 2):
            # print('j', j)
            # fasta indexing starts at 1
            if aa[j + 1] != sequence[j + 1 + shift]:
                all_match = False 
                break 
        if all_match:
            print('shift is:', shift)
            return shift
            break
    if not all_match:
        return 0
        print('no matching shift found.')
        raise Exception ('no matching shift found.')

# parse pdb file
parser = PDBParser()
model = parser.get_structure("", PDB_FILE)[0]

dssp = DSSP(model, PDB_FILE)

# DSSP data key is a tuple (chain_id, (?, res_id, ?))
# DSSP data value is DSSP index, aa, sec struct, rsa, ...
# rsa_positions = [key[1][1] for key in list(dssp.keys())]
rsa_value = [dssp[key][3] for key in list(dssp.keys())]
rsa_aa = [dssp[key][1] for key in list(dssp.keys())]
# PDB indexing!!! check if indexes work with fasta indices, check with wtSequence._data

# if I don't catch it will show up as exception in dashboard 
# try:
shift = check_for_shift(wtSequence, rsa_aa)
# except Exception as e:
#     print('Exception for matching pdb sites with fasta sites:', e)

rsa_positions = list(range(1 + shift, len(rsa_aa) + 1 + shift))

def make_fig_cons_rsa():
    print("making fig_cons_rsa")
    fig_cons_rsa = plotly.tools.make_subplots(
        shared_xaxes=True, vertical_spacing=0.08, 
        rows=4, cols=1, row_heights=[0.15, 0.15, 0.3, 0.4], subplot_titles=("(A) Coservation of homologous sequence", "(B) Relative solvent accessibility", "(C) {} of variants with mutation at respective site".format(ATTRIBUTE_COL), "(D) Median {} of variants with respective mutation at respective site".format(ATTRIBUTE_COL)))
    # fig.add_trace(px.bar(conservation_df, x='site' , y='conservation', color='gaps'), row=1, col=1)
    fig_cons_rsa.append_trace(
        go.Bar(
            x = conservation_df['site'], y=conservation_df['conservation'], marker=dict(color=conservation_df['gaps'], colorbar=dict(title='% gaps', y=1.02, len=0.15, yanchor="top"))
        ), row=1,col=1)
    fig_cons_rsa.update_yaxes(row=1, col=1, title_text = "information content")

    fig_cons_rsa.append_trace(
        go.Bar(x=rsa_positions, y=rsa_value), row=2, col=1)

    # same limits as conservation plot, over whole sequence
    fig_cons_rsa.update_xaxes(range=[1, seq_length + 0.5], row=2, col=1)
    fig_cons_rsa.update_yaxes(row=2, col=1, title_text = "rsa")
    
    fig_cons_rsa.update_layout(showlegend=False, height=800, margin=dict(b=20, t=30))
    return fig_cons_rsa



def apply_filter(df, attribute, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa):
    """apply filter to dataframe, possibillity to apply only selected filters by passing in None"""

    df_filtered = df.copy()

    if numMutation is not None:
        df_filtered = df_filtered[df_filtered["numMutation"] >= numMutation]

    if attribute is not None:
        df_filtered = df_filtered[df_filtered[ATTRIBUTE_COL] >= attribute[0]]
        df_filtered = df_filtered[df_filtered[ATTRIBUTE_COL] <= attribute[1]]

    if sites_mutation is not None:
        if sites_mutation != []:
            matching_variants = []
            for site in sites_mutation:
                for i, row in df_filtered.iterrows():
                    if site in [int(mutation[1:-1]) for mutation in row["mutationAll"]]:
                        matching_variants.append(row[VARIANT_COL])
                    elif site == 'None' and row["mutationAll"] == []:
                        # print('site is None', row["mutationAll"])
                        matching_variants.append(row[VARIANT_COL])

            matching_variants = list(set(matching_variants))

            df_filtered = df_filtered[df_filtered[VARIANT_COL].isin(matching_variants)]

    if sites_new_mutation is not None:
        if sites_new_mutation != []:
            matching_variants = []
            for site in sites_new_mutation:
                for i, row in df_filtered.iterrows():
                    # check if variant has new mutation at all
                    if site == 'None' and row[MUTATION_COL] == None:
                        # print('new_mutation is none', row[MUTATION_COL] == None)
                        matching_variants.append(row[VARIANT_COL])
                    elif type(row[MUTATION_COL]) == str and site == int(
                        row[MUTATION_COL][1:-1]
                    ):
                        matching_variants.append(row[VARIANT_COL])

            matching_variants = list(set(matching_variants))

            df_filtered = df_filtered[df_filtered[VARIANT_COL].isin(matching_variants)]

    if children is not None:
        # filter every variant that has that variant Id or has it in it's parentAll
        # if None selected: don't filter
        if children != 'None':
            # print('children', children)
            matching_variants = []
            for i, row in df_filtered.iterrows():
                # print('row parentAll', row['parentAll'])
                if children in row["parentAll"]:
                    # print('is in row')
                    matching_variants.append(row[VARIANT_COL])

            matching_variants = list(set(matching_variants))

            df_filtered = df_filtered[df_filtered[VARIANT_COL].isin(matching_variants)]
    
    if exact_site != None and exact_aa != None:
        matching_variants = []
        for i, row in df_filtered.iterrows():
            # print('row parentAll', row['parentAll'])
            if row["sequence"][exact_site - 1] == exact_aa:
                # print('is in row')
                matching_variants.append(row[VARIANT_COL])

        matching_variants = list(set(matching_variants))

        df_filtered = df_filtered[df_filtered[VARIANT_COL].isin(matching_variants)]

    return df_filtered



def populate_mutation_site_dropdown(df, return_type):
    """populate options and values of mutation sites dropdown"""

    mutation_list = []
    for row in [row['mutationAll'] for i, row in df.iterrows() if row["mutationAll"] != []]:
        for mutation in row:
            mutation_list.append(mutation)

    site_list = sorted(list(set([int(mutation[1:-1]) for mutation in mutation_list])))
    site_list.append('None')

    if return_type == 'options':
        return [
        {"label": site, "value": site}
        for site in site_list]

    if return_type == 'value':
        return site_list



def populate_new_mutation_site_dropdown(df, return_type):
    """populate options and values of new mutation sites dropdown"""

    site_list = sorted(list(set([int(row[MUTATION_COL][1:-1]) for i, row in df.iterrows() if row[MUTATION_COL] != None])))
    site_list.append('None')

    if return_type == 'options':
        return [
        {"label": site, "value": site}
        for site in site_list]

    if return_type == 'value':
        return site_list


def populate_children_options(df):
    '''returns label value dict with all variants that have children and additionally one None option '''

    all_variants = df['variant'].tolist()
    variants_with_children = []

    for variant in all_variants:
        for i, row in df.iterrows():
            if variant in row['parentAll']:
                # this will prevent duplications while keeping order in place
                if variant not in variants_with_children:
                    variants_with_children.append(variant)
                    continue

    return [{"label": variant, "value": variant} for variant in variants_with_children]# + [{"label": 'None', "value": "None"}]



#########################
# Dashboard Layout / View
#########################

# columns I want to display in the table
columns_display = {
    VARIANT_COL: "text",
    PARENT_COL: "text",
    MUTATION_COL: "text",
    ATTRIBUTE_COL: "numeric",
    "numMutation": "numeric",
}

# Set up Dashboard and create layout
app = dash.Dash(
    __name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"]
)
# app = dash.Dash(__name__)

app.config["suppress_callback_exceptions"] = True


app.layout = html.Div(
    [
        # Page header
        html.Div([html.H1("Protein Engineering Dashboard")]),
        html.H3("Filter Data"),
        # filter options
        html.Div(
            [
                # multi dropdown filter options
                html.Label("select filter options"),
                dcc.Dropdown(
                    id="filter-option-dropdown",
                    # populatte options with all available sites
                    options=[
                        {"label": ATTRIBUTE_COL, "value": "attribute-range-slider"},
                        {"label": "n mutations", "value": "numMutation-input"},
                        {"label": "mutation sites", "value": "mutation-site-dropdown"},
                        {
                            "label": "new mutation sites",
                            "value": "new-mutation-site-dropdown",
                        },
                        {"label": "children", "value": "children-dropdown"},
                        {"label": "exact position", "value": "exact-position"},
                    ],
                    placeholder="select filter",
                    multi=True,
                ),
            ]
        ),
        html.Div(
            [
                html.Div(
                    [
                        # range slider for attribute
                        html.Label("Select {} range".format(ATTRIBUTE_COL)),
                        html.Div(id="attribute-slider-histogram-container"),
                        dcc.RangeSlider(
                            id="attribute-range-slider",
                            min=df_extended[ATTRIBUTE_COL].min(),
                            max=df_extended[ATTRIBUTE_COL].max(),
                            step=0.1,
                            value=[df_extended[ATTRIBUTE_COL].min(), df_extended[ATTRIBUTE_COL].max()],
                            marks={
                                df_extended[ATTRIBUTE_COL].min(): str(df_extended[ATTRIBUTE_COL].min()),
                                df_extended[ATTRIBUTE_COL].max(): str(df_extended[ATTRIBUTE_COL].max()),
                            },
                        ),
                    ],
                    id="attribute-range-slider-container",
                    style={
                        "display": "none",
                        "padding-bottom": "30px",
                        "break-inside": "avoid-column",
                        "column-break-inside": "avoid",
                    },
                ),  # ????
                html.Div(
                    [
                        # number input numMutation
                        html.Label("Select minimum number of mutations"),
                        dcc.Input(
                            id="numMutation-input",
                            type="number",
                            # debounce=True,  # limits firing after every change, does not change until e.g. hit enter
                            min=df_extended["numMutation"].min(),
                            max=df_extended["numMutation"].max(),
                            step=1,
                            value=df_extended["numMutation"].min(),
                        ),
                    ],
                    id="numMutation-input-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    [
                        # multi dropdown mutation site
                        html.Label(
                            "select mutation site(s) (in relation to wt sequence)"
                        ),
                        dcc.Dropdown(
                            id="mutation-site-dropdown",
                            # populatte options with all available sites
                            options=populate_mutation_site_dropdown(df_extended, 'options'),
                            # value=populate_mutation_site_dropdown(df_extended, 'value'),
                            placeholder='Select mutation',
                            multi=True,
                        ),
                    ],
                    id="mutation-site-dropdown-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    [
                        # multi dropdown mutation site
                        html.Label("select site of new mutation"),
                        dcc.Dropdown(
                            id="new-mutation-site-dropdown",
                            # populatte options with all available sites, only filtering changes
                            options=populate_new_mutation_site_dropdown(df_extended, 'options'),
                            # value=populate_new_mutation_site_dropdown(df_extended, 'value'),
                            placeholder='Select mutation',
                           multi=True,
                        ),
                    ],
                    id="new-mutation-site-dropdown-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    [
                        # dropdown children of variant
                        html.Label(
                            "select variant whose children to select"
                        ),
                        dcc.Dropdown(
                            id="children-dropdown",
                            # populate options with all available variants
                            options=populate_children_options(df_extended),
                            # pick wt for initial value
                            # value='None',
                            # multi=False,
                        ),
                    ],
                    id="children-dropdown-container",
                    style={"display": "none"},
                ),
                html.Div(
                    [
                        # 2 dropdowns for exact position
                        html.Label(
                            "Select position and amino acid"
                        ),
                        dcc.Dropdown(
                            id="exact-position-site",
                            # populate options with all available sites
                            options=[{"label": str(site), "value": site} for site in list(range(1, seq_length + 1))],
                            placeholder='Select position', 
                        ),
                        dcc.Dropdown(
                            id="exact-position-aa",
                            # populate options with available aa at that site with callback, none if position is none
                            placeholder='Select amino acid', 
                        ),
                    ],
                    id="exact-position-container",
                    style={"display": "none"},
                ),
                html.Div(id="filtered-data", style={"display": "none"}),
            ],
            style={"padding-bottom": "30px", "columnCount": 1},
        ),
        # html.Div(id="filter-values-container", style={"padding-bottom": "10px"}),
        # one container for table and graph so I can fit them next to each other
        html.Div(
            [
                # container for table
                html.Div(id="filtered-data-table-container", className="five columns"),
                # container for plot
                html.Div(id="variant-plot-container", className="seven columns"),
                # paddinf as large as table
            ],
            style={"padding-bottom": "520px"},
        ),
        html.Br(),
        html.Div(
            [
                dcc.RadioItems(
                    id="violin-display-option",
                    options=[
                        {"label": "box-plot", "value": "box-plot"},
                        {"label": "violin", "value": "violin"},
                        # {'label': 'bar', 'value': 'bar'},
                        {"label": "scatter", "value": "scatter"},
                    ],
                    value="box-plot",
                    labelStyle={"display": "inline-block"},
                ),
            ],
            style={"text-align": "center"},
        ),
        html.Div(id="violin-plot-container", className="twelve-columns"),
        # html.Hr(),
        # html.Hr(),
        html.H3(["Define Hit"], style={"padding-top": "20px"}),
        # hit definition filter
        html.Div(
            [
                # multi dropdown filter options
                html.Label("select hit definition options"),
                dcc.Dropdown(
                    id="hit-option-dropdown",
                    # populate options with all available sites
                    options=[
                        {"label": ATTRIBUTE_COL, "value": "attribute-range-slider-hit"},
                        {"label": "n mutations", "value": "numMutation-input-hit"},
                        {
                            "label": "mutation sites",
                            "value": "mutation-site-dropdown-hit",
                        },
                        {
                            "label": "new mutation sites",
                            "value": "new-mutation-site-dropdown-hit",
                        },
                        {"label": "children", "value": "children-dropdown-hit"},
                        {"label": "exact position", "value": "exact-position-hit"},
                    ],
                    placeholder="select filter",
                    multi=True,
                ),
            ]
        ),
        html.Div(
            [
                html.Div(
                    # range slider for attribute
                    id="attribute-range-slider-hit-container",
                    style={"display": "none", "padding-bottom": "30px"},
                ),  # ????
                html.Div(
                    # number input numMutation
                    id="numMutation-input-hit-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    # multi dropdown mutation site
                    id="mutation-site-dropdown-hit-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    # multi dropdown mutation site
                    id="new-mutation-site-dropdown-hit-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    # dropdown children
                    id="children-dropdown-hit-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(
                    # drop down exact position
                    id="exact-position-hit-container",
                    style={"display": "none", "break-after": "column"},
                ),
                html.Div(id="hit-data", style={"display": "none"}),
            ],
            style={"padding-bottom": "30px", "columnCount": 1},
        ),
        # html.Div(id="hit-values-container", style={"padding-bottom": "10px"}),
        html.H5("Sequence logo of homologous sequences (HHblits, ggseqlogo)"),
        html.Div(html.Img(src=calculate_homology_logo()),
            style={"overflow": "auto", "paddingRight": "20px", "padding-left": "20px"},
        ),

        html.Div(id="seq-logos-container"),
        html.Div(id="swarm-container"),
    ],
    style={"columnCount": 1},
)


#############################################
# Interaction Between Components / Controller
#############################################

################################# Filtering ####################################################

@app.callback(
    [
        Output("attribute-range-slider-container", "style"),
        Output("attribute-range-slider", "value"),
        Output("numMutation-input-container", "style"),
        Output("numMutation-input", "value"),
        Output("mutation-site-dropdown-container", "style"),
        Output("mutation-site-dropdown", "value"),
        Output("new-mutation-site-dropdown-container", "style"),
        Output("new-mutation-site-dropdown", "value"),
        Output("children-dropdown-container", "style"),
        Output("children-dropdown", "value"),
        Output("exact-position-container", "style"),
        Output("exact-position-site", "value"),
        Output("exact-position-aa", "value"),
    ],
    [Input("filter-option-dropdown", "value"),],
    # prevent_initial_call=True
)
def get_selected_filter(value):
    """handle visibility of filter options, reset values if unselected"""

    # if value is not None:
    return_list = []

    # must remain in correct order like output, maybe list is safer
    filter_options = [
        "attribute-range-slider",
        "numMutation-input",
        "mutation-site-dropdown",
        "new-mutation-site-dropdown",
        "children-dropdown"
    ]
    filter_default = {
        "attribute-range-slider": [df_extended[ATTRIBUTE_COL].min(), df_extended[ATTRIBUTE_COL].max()],
        "numMutation-input": df_extended["numMutation"].min(),
        "mutation-site-dropdown": None,
        "new-mutation-site-dropdown": None,
        "children-dropdown": None
    }

    for filter_option in filter_default:
        if value is None or filter_option not in value:
            return_list.append({"display": "none"})
            return_list.append(filter_default[filter_option])
        else:
            return_list.append({"display": "initial"})
            return_list.append(dash.no_update)

    if value is None or 'exact-position' not in value:
        return_list.append({"display": "none"})
        return_list.append(None)
        return_list.append(None)
    else:
        return_list.append({"display": "initial"})
        return_list.append(dash.no_update)
        return_list.append(dash.no_update)
 
    return return_list


@app.callback(
    Output('exact-position-aa', 'options'),
    [Input('exact-position-site', 'value')]
)
def populate_exact_position_aa(site):
    '''return all available aa at selected site'''
    # print('type site', type(site))
    if site is not None:
        available_aa = list(set([sequence[site - 1] for sequence in df_extended['sequence']]))
        dict_to_return = [{"label": aa, "value": aa} for aa in available_aa]
        return dict_to_return

    raise dash.exceptions.PreventUpdate



@app.callback(
    Output("attribute-slider-histogram-container", "children"),
    [
        Input("numMutation-input", "value"),
        Input("mutation-site-dropdown", "value"),
        Input("new-mutation-site-dropdown", "value"),
        Input("children-dropdown", "value"),
        Input("exact-position-site", "value"),
        Input("exact-position-aa", "value"),
    ],
)
def draw_attribute_hist(numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa):
    """histogram on top of attribute slider"""

    df_filtered = apply_filter(
        df_extended, None, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa
    )

    fig = px.histogram(
        df_filtered,
        x=ATTRIBUTE_COL,
        nbins=100,
        template="plotly_white",
        color_discrete_sequence=["gray"],
    )
    fig.update_xaxes(
        showticklabels=True,
        color="#fff",
        title_text="",
        gridcolor="#fff",
        range=[df_extended[ATTRIBUTE_COL].min(), df_extended[ATTRIBUTE_COL].max()],
    )
    fig.update_yaxes(
        showticklabels=False, color="#fff", title_text="", gridcolor="#fff"
    )
    fig.update_layout(margin={"b": 5, "l": 0, "r": 23, "t": 0, "pad": 0}, height=50)

    return dcc.Graph(id="attribute-slider-histogram", figure=fig)


@app.callback(
    Output("filtered-data", "children"),
    [
        Input("attribute-range-slider", "value"),
        Input("numMutation-input", "value"),
        Input("mutation-site-dropdown", "value"),
        Input("new-mutation-site-dropdown", "value"),
        Input("children-dropdown", "value"),
        Input("exact-position-site", "value"),
        Input("exact-position-aa", "value"),
    ],
)
def filter_data(attribute, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa):
    """apply filter, store json in hidden Div"""

    df_filtered = apply_filter(
        df_extended, attribute, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa
    )

    return df_filtered.to_json(date_format="iso", orient="split")




@app.callback(
    Output("filtered-data-table-container", "children"),
    [Input("filtered-data", "children"),],
)
def filter_data(filtered_data):
    """render sortable table with filtered data"""

    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)

    return dash_table.DataTable(
        id="data-overview-table",
        columns=[
            {"name": i, "id": i, "type": columns_display[i]} for i in columns_display
        ],
        # data will be updated by callback
        data=df_filtered.loc[:, columns_display.keys()].to_dict("records"),
        # filter_action='native',
        sort_action="native",
        style_table={"max-height": "500px", "overflowY": "auto"},
    )


################################ Hit Definition ################################################




@app.callback(
    [
        Output("attribute-range-slider-hit-container", "children"),
        Output("numMutation-input-hit-container", "children"),
        Output("mutation-site-dropdown-hit-container", "children"),
        Output("new-mutation-site-dropdown-hit-container", "children"),
        Output("children-dropdown-hit-container", "children"),
        Output("exact-position-hit-container", "children"),
    ],
    [Input("filtered-data", "children"),],
)
def populate_hit_filter(filtered_data):
    """populate hit definition filter"""

    print("** populating hit options")

    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)

    # print('df_filtered', df_filtered.iloc[0])

    populated_hit_filter = []
    return [
        html.Div(
            [
                html.Label("Select {} range".format(ATTRIBUTE_COL)),
                html.Div(id="attribute-slider-histogram-hit-container"),
                dcc.RangeSlider(
                    id="attribute-range-slider-hit",
                    min=df_filtered[ATTRIBUTE_COL].min(),
                    max=df_filtered[ATTRIBUTE_COL].max(),
                    step=0.1,
                    value=[df_filtered[ATTRIBUTE_COL].min(), df_filtered[ATTRIBUTE_COL].max()],
                    marks={
                        df_filtered[ATTRIBUTE_COL].min(): str(df_filtered[ATTRIBUTE_COL].min()),
                        df_filtered[ATTRIBUTE_COL].max(): str(df_filtered[ATTRIBUTE_COL].max()),
                    },
                ),
            ]
        ),
        html.Div(
            [
                html.Label("Select minimum number of mutations"),
                dcc.Input(
                    id="numMutation-input-hit",
                    type="number",
                    # debounce=True,  # limits firing after every change, does not change until e.g. hit enter
                    min=df_filtered["numMutation"].min(),
                    max=df_filtered["numMutation"].max(),
                    step=1,
                    value=df_filtered["numMutation"].min(),
                ),
            ]
        ),
        html.Div(
            [
                html.Label("select mutation site(s)"),
                dcc.Dropdown(
                    id="mutation-site-dropdown-hit",
                    # populatte options with all available sites
                    options=populate_mutation_site_dropdown(df_filtered, 'options'),
                    # value=populate_mutation_site_dropdown(df_filtered, 'value'),
                    multi=True,
                ),
            ]
        ),
        html.Div(
            [
                html.Label("select mutation site(s)"),
                dcc.Dropdown(
                    id="new-mutation-site-dropdown-hit",
                    # populatte options with all available sites
                    options=populate_new_mutation_site_dropdown(df_filtered, 'options'),
                    # value=populate_new_mutation_site_dropdown(df_filtered, 'value'),
                    multi=True,
                ),
            ]
        ),
        html.Div(
            [
                html.Label("select variant whose children to select"),
                dcc.Dropdown(
                    id="children-dropdown-hit",
                    # populatte options with all available sites
                    options=populate_children_options(df_filtered),
                    # value='None',
                ),
            ]
        ),
        html.Div(
            [
                # 2 dropdowns for exact position
                html.Label(
                    "Select position and amino acid"
                ),
                dcc.Dropdown(
                    id="exact-position-site-hit",
                    # populate options with all available sites
                    options=[{"label": str(site), "value": site} for site in list(range(1, seq_length + 1))],
                    placeholder='Select position', 
                ),
                dcc.Dropdown(
                    id="exact-position-aa-hit",
                    # populate options with available aa at that site with callback, none if position is none
                    placeholder='Select amino acid', 
                ),
            ],
        ),

    ]




@app.callback(
    [
        Output("attribute-range-slider-hit-container", "style"),
        Output("attribute-range-slider-hit", "value"),
        Output("numMutation-input-hit-container", "style"),
        Output("numMutation-input-hit", "value"),
        Output("mutation-site-dropdown-hit-container", "style"),
        Output("mutation-site-dropdown-hit", "value"),
        Output("new-mutation-site-dropdown-hit-container", "style"),
        Output("new-mutation-site-dropdown-hit", "value"),
        Output("children-dropdown-hit-container", "style"),
        Output("children-dropdown-hit", "value"),
        Output("exact-position-hit-container", "style"),
        Output("exact-position-site-hit", "value"),
        Output("exact-position-aa-hit", "value"),
    ],
    [Input("hit-option-dropdown", "value"), Input("filtered-data", "children"),],
    prevent_initial_call=True,
)
def get_selected_filter(value, filtered_data):
    """handle visibility of hit options, reset values if unselected"""

    # print("visibility hit filter")
    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)

    # if value is not None:
    return_list = []

    # must remain in correct order like output, maybe list is safer
    filter_options = [
        "attribute-range-slider-hit",
        "numMutation-input-hit",
        "mutation-site-dropdown-hit",
        "new-mutation-site-dropdown-hit",
        "children-dropdown-hit",
    ]
    filter_default = {
        "attribute-range-slider-hit": [df_filtered[ATTRIBUTE_COL].min(), df_filtered[ATTRIBUTE_COL].max()],
        "numMutation-input-hit": df_filtered["numMutation"].min(),
        "mutation-site-dropdown-hit":  None,
        "new-mutation-site-dropdown-hit":  None,
        "children-dropdown-hit": None
    }

    for filter_option in filter_default:
        if value is None or filter_option not in value:
            return_list.append({"display": "none"})
            return_list.append(filter_default[filter_option])
        else:
            return_list.append({"display": "initial"})
            return_list.append(dash.no_update)

    if value is None or 'exact-position-hit' not in value:
        return_list.append({"display": "none"})
        return_list.append(None)
        return_list.append(None)
    else:
        return_list.append({"display": "initial"})
        return_list.append(dash.no_update)
        return_list.append(dash.no_update)

    return return_list


@app.callback(
    Output('exact-position-aa-hit', 'options'),
    [Input('exact-position-site-hit', 'value'), 
    Input("filtered-data", "children")]
)
def populate_exact_position_aa_hit(site, filtered_data):
    '''return all available aa at selected site'''
    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)

    # print('type site', type(site))
    if site is not None:
        available_aa = list(set([sequence[site - 1] for sequence in df_filtered['sequence']]))
        # print('available aa:', available_aa)
        dict_to_return = [{"label": aa, "value": aa} for aa in available_aa]
        # print('dict to return', dict_to_return)
        return dict_to_return
    raise dash.exceptions.PreventUpdate


@app.callback(
    Output("attribute-slider-histogram-hit-container", "children"),
    [
        Input("numMutation-input-hit", "value"),
        Input("mutation-site-dropdown-hit", "value"),
        Input("new-mutation-site-dropdown-hit", "value"),
        Input("children-dropdown-hit", "value"),
        Input("exact-position-site-hit", "value"),
        Input("exact-position-aa-hit", "value"),
    ],
)
def draw_attribute_hist_hit(numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa):
    """histogram on top of attribute slider"""

    df_hit = apply_filter(
        df_extended, None, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa
    )

    fig = px.histogram(
        df_hit,
        x=ATTRIBUTE_COL,
        nbins=100,
        template="plotly_white",
        color_discrete_sequence=["gray"],
    )
    fig.update_xaxes(
        showticklabels=True,
        color="#fff",
        title_text="",
        gridcolor="#fff",
        range=[df_extended[ATTRIBUTE_COL].min(), df_extended[ATTRIBUTE_COL].max()],
    )
    fig.update_yaxes(
        showticklabels=False, color="#fff", title_text="", gridcolor="#fff"
    )
    fig.update_layout(margin={"b": 5, "l": 0, "r": 23, "t": 0, "pad": 0}, height=50)

    return dcc.Graph(id="attribute-slider-histogram-hit", figure=fig)


@app.callback(
    Output("hit-data", "children"),
    [
        Input("attribute-range-slider-hit", "value"),
        Input("numMutation-input-hit", "value"),
        Input("mutation-site-dropdown-hit", "value"),
        Input("new-mutation-site-dropdown-hit", "value"),
        Input("children-dropdown-hit", "value"),
        Input("exact-position-site-hit", "value"),
        Input("exact-position-aa-hit", "value"),
    ],
)
def hit_data(attribute, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa):
    """apply filter, store json in hidden Div"""

    # print(ATTRIBUTE_COL, attribute, 'numMutation', numMutation, 'sites_mutation', sites_mutation, 'sites_new_mutation', sites_new_mutation)
    df_hit = apply_filter(
        df_extended, attribute, numMutation, sites_mutation, sites_new_mutation, children, exact_site, exact_aa
    )
    # print('df_hit before saving to json', df_hit)

    return df_hit.to_json(date_format="iso", orient="split")


################################# plots ################################################################

@app.callback(
    Output("variant-plot-container", "children"), [Input("filtered-data", "children"),]
)
def update_graphs(filtered_data):
    """Variant plot, updates according to filter selection in table"""
    # https://dash.plotly.com/datatable/interactivity

    dff = pd.read_json(filtered_data, orient="split")
    dff = dff.where(pd.notnull(dff), None)

    filter_applied = True if len(df_extended) != len(dff) else False

    fig = go.Figure()

    # draw line from each variant to it's parents, shows opaque lines for all variants even if filtered

    # only draw transparent lines if something was filtered
    if filter_applied:
        # one line/trace for each variant and each parent
        for i, row in df_extended.iterrows():
            if not pd.isnull(row[PARENT_COL]):
                for parent in row[PARENT_COL].split(" "):
                    x = [
                        row["numMutation"],
                        df_extended.loc[df_extended[VARIANT_COL] == parent][
                            "numMutation"
                        ].tolist()[0],
                    ]
                    y = [
                        row[ATTRIBUTE_COL],
                        df_extended.loc[df_extended[VARIANT_COL] == parent][
                            ATTRIBUTE_COL
                        ].tolist()[0],
                    ]
                    fig.add_trace(
                        go.Scatter(
                            x=x,
                            y=y,
                            mode="lines",
                            line=dict(width=1, color="rgba(30, 50, 180, 0.2)"),
                            # opacity=0.2,
                            name="",
                        )
                    )

    # all variants that were selected
    parentVariantsToExplore = set(dff[VARIANT_COL].tolist())
    parentVariantsExplored = set()

    # trace back path for each variant that is selected
    # as long as there are parents to explore
    while parentVariantsToExplore != set():
        # remove one item
        parentVariant = parentVariantsToExplore.pop()
        # print('variant', parentVariant)
        i = df_extended[df_extended[VARIANT_COL] == parentVariant].index.tolist()[0]

        parentVariantsExplored.add(parentVariant)

        if not pd.isnull(df_extended[PARENT_COL].iloc[i]):
            for parent in df_extended[PARENT_COL].iloc[i].split(" "):

                x = [
                    df_extended["numMutation"].iloc[i],
                    df_extended.loc[df_extended[VARIANT_COL] == parent][
                        "numMutation"
                    ].tolist()[0],
                ]
                y = [
                    df_extended[ATTRIBUTE_COL].iloc[i],
                    df_extended.loc[df_extended[VARIANT_COL] == parent][ATTRIBUTE_COL].tolist()[
                        0
                    ],
                ]
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(width=1.5, color="rgba(30, 50, 180, 0.5)"),
                        # opacity=0.5,
                        name="",
                    )
                )

                if parent not in parentVariantsExplored:
                    parentVariantsToExplore.add(parent)

    # shows points for selected variants and all variants with different opacity
    # one trace for selected variants and one for all
    fig.add_trace(
        go.Scatter(
            x=dff["numMutation"],
            y=dff[ATTRIBUTE_COL],
            mode="markers",
            opacity=1,
            marker_color="rgba(30, 50, 180, 1)",
            name="",
        )
    )

    text = [
        "Variant: {}<br>Mutation: {}".format(row[VARIANT_COL], row[MUTATION_COL])
        for i, row in df_extended.iterrows()
    ]

    fig.add_trace(
        go.Scatter(
            x=df_extended["numMutation"],
            y=df_extended[ATTRIBUTE_COL],
            mode="markers",
            # opacity=0.2,
            marker_color="rgba(30, 50, 180, 0.2)",
            text=text,
            hovertemplate="%{text}<br>n mutations: %{x}<br>" + ATTRIBUTE_COL + ": %{y:.1f}",
            name="",
        )
    )

    fig.update_layout(
        showlegend=False,
        title=dict(text="Variant Plot", x=0.5, y=0.9, xanchor="center", yanchor="top"),
        xaxis_title="Number of Mutations",
        yaxis_title=ATTRIBUTE_COL,
        font=dict(size=14),
    )

    return dcc.Graph(id="variant-plot", figure=fig)



@app.callback(
    Output("violin-plot-container", "children"),
    [Input("filtered-data", "children"), Input("violin-display-option", "value")],
)
def update_violinplot(filtered_data, display_option):
    """violin plot, recieves filtered data as json"""

    fig_cons_rsa = make_fig_cons_rsa()

    dff = pd.read_json(filtered_data, orient="split")
    dff = dff.where(pd.notnull(dff), None)

    df_violin = pd.DataFrame(columns=[VARIANT_COL, "site", ATTRIBUTE_COL])

    sites = list(range(1, len(dff["sequence"].iloc[0]) + 1))

    for i, row in dff.iterrows():
        for mutation in row["mutationAll"]:
            df_violin = df_violin.append(
                {
                    VARIANT_COL: row[VARIANT_COL],
                    "site": int(mutation[1:-1]),
                    ATTRIBUTE_COL: row[ATTRIBUTE_COL],
                },
                ignore_index=True,
            )

    if display_option == "violin":

        for site in sites:
            fig_cons_rsa.append_trace(
                go.Violin(
                    x=df_violin["site"][df_violin["site"] == site],
                    y=df_violin[ATTRIBUTE_COL][df_violin["site"] == site],
                    text=df_violin[VARIANT_COL][df_violin["site"] == site],
                    hovertemplate="Variant: %{text}<br>"
                    + ATTRIBUTE_COL + ": %{y}<br>"
                    + "Position: %{x}",
                    name="",  # trace name
                    box_visible=False,
                    meanline_visible=True,
                    # points='all',
                    pointpos=0,
                    jitter=0.5,
                    line_color="rgba(30, 50, 180, 1)",
                    line_width=0,
                    # scalemode='count',
                    spanmode="hard",
                    width=0.8,
                    # opacity=0.3
                )
                , row=3, col=1
            )

    elif display_option == "box-plot":
        # fig = px.box(df_violin, x="site", y=ATTRIBUTE_COL)
        # fig.update_layout(margin={"t": 100})
        # fig.append_trace(px.box(df_violin, x="site", y=ATTRIBUTE_COL), row=3, col=1)
        for site in sites:
            fig_cons_rsa.append_trace(
                go.Box(
                    y=df_violin[ATTRIBUTE_COL][df_violin["site"] == site],
                    x=df_violin["site"][df_violin["site"] == site],
                    marker={'color': "rgba(30, 50, 180, 1)"},
                )
                , row=3, col=1
            )

        # fig_cons_rsa.update_layout(margin={"t": 100})

    elif display_option == "scatter":
        # fig = go.Figure()
        fig_cons_rsa.append_trace(
            go.Scatter(
                x=df_violin["site"],
                y=df_violin[ATTRIBUTE_COL],
                mode="markers",
                marker=dict(
                    color="rgba(30, 50, 180, 0.5)",
                    size=5
                ),
                text=df_violin[VARIANT_COL],
                hovertemplate="variant: %{text}<br>"
                + ATTRIBUTE_COL + ": %{y}<br>"
                + "Position: %{x}",
                name="",
            )
        , row=3, col=1)

    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)

    len_seq = len(wtSequence)
    sites = list(range(1, len(df_filtered["sequence"].iloc[0]) + 1))

    matrix = [
        [
            []
            for site in sites
        ]
        for aa in aa_alphabet
    ]

    for i, row in df_filtered.iterrows():
        if row["mutationAll"] == None:
            continue 
        for mutation in row["mutationAll"]:
            matrix[aa_alphabet.index(mutation[-1])][int(mutation[1:-1])].append(row[ATTRIBUTE_COL])

    for aa in range(len(aa_alphabet)):
        for i in range(len(sites)):
            matrix[aa][i] = statistics.median(matrix[aa][i]) if matrix[aa][i] != [] else None

    fig_cons_rsa.append_trace(go.Heatmap(
        z=matrix,
        x=sites,
        y=aa_alphabet,
        colorbar=dict(title=ATTRIBUTE_COL, y=-0.02, len=0.35, yanchor="bottom"),
        hovertemplate="site: %{x}<br>new aa: %{y}<br>" + ATTRIBUTE_COL + ": %{z:.1f}",
    ), row=4, col=1)

    fig_cons_rsa.update_xaxes(
        showgrid=True,
        range=[list(sites)[0], list(sites)[-1]],
        tickvals=list(range(5, list(sites)[-1] + 1, 5)),
        row=3, col=1
    )
    fig_cons_rsa.update_yaxes(range=[df_violin[ATTRIBUTE_COL].min() - 1, df_violin[ATTRIBUTE_COL].max() + 1], title_text = ATTRIBUTE_COL, row=3, col=1)

    fig_cons_rsa.update_xaxes(
        title_text = "position",
        tickvals=list(range(5, list(sites)[-1] + 1, 5)),
        row=4, col=1
    )
    fig_cons_rsa.update_yaxes(
        tickmode="linear",
        title_text="amino acid",
        row=4, col=1
    )

    return dcc.Graph(id="violin-plot", figure=fig_cons_rsa)


@app.callback(
    Output("swarm-container", "children"),
    [Input("filtered-data", "children"), Input("hit-data", "children"),],
)
def update_swarm(filtered_data, hit_data):
    if hit_data is None:
        print("not updating swarm")
        raise dash.exceptions.PreventUpdate

    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)
    df_hit = pd.read_json(hit_data, orient="split")
    df_hit = df_hit.where(pd.notnull(df_hit), None)
    df_not_hit = df_filtered[df_filtered['variant'].isin(set(df_filtered['variant']) - set(df_hit['variant']))]

    if len(df_filtered) == len(df_hit):
        raise dash.exceptions.PreventUpdate

    df_swarm_hit = pd.DataFrame(columns=[VARIANT_COL, "site", ATTRIBUTE_COL])
    df_swarm_non_hit = pd.DataFrame(columns=[VARIANT_COL, "site", ATTRIBUTE_COL])

    sites = list(range(1, len(df_filtered["sequence"].iloc[0]) + 1))

    for i, row in df_hit.iterrows():
        for mutation in row["mutationAll"]:
            df_swarm_hit = df_swarm_hit.append(
                {
                    VARIANT_COL: row[VARIANT_COL],
                    "site": int(mutation[1:-1]),
                    ATTRIBUTE_COL: row[ATTRIBUTE_COL],
                },
                ignore_index=True,
            )
    for i, row in df_not_hit.iterrows():
        for mutation in row["mutationAll"]:
            df_swarm_non_hit = df_swarm_non_hit.append(
                {
                    VARIANT_COL: row[VARIANT_COL],
                    "site": int(mutation[1:-1]),
                    ATTRIBUTE_COL: row[ATTRIBUTE_COL],
                },
                ignore_index=True,
            )


    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df_swarm_hit["site"],
            y=df_swarm_hit[ATTRIBUTE_COL],
            mode="markers",
            # marker_color="blue",
            marker=dict(
                color="rgba(214, 39, 40, 0.5)",
                size=5
            ),
            # opacity=0.3,
            text=df_swarm_hit[VARIANT_COL],
            hovertemplate="variant: %{text}<br>"
            + ATTRIBUTE_COL + ": %{y}<br>"
            + "Position: %{x}",
            name="hit",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=df_swarm_non_hit["site"],
            y=df_swarm_non_hit[ATTRIBUTE_COL],
            mode="markers",
            # marker_color="blue",
            marker=dict(
                color="rgba(50, 30, 180, 0.5)",
                size=5
            ),
            # opacity=0.3,
            text=df_swarm_non_hit[VARIANT_COL],
            hovertemplate="variant: %{text}<br>"
            + ATTRIBUTE_COL + ": %{y}<br>"
            + "Position: %{x}",
            name="non-hit",
        )
    )

    fig.update_xaxes(
        title_text = "position",
        showgrid=True,
        range=[list(sites)[0], list(sites)[-1]],
        tickvals=list(range(5, list(sites)[-1] + 1, 5)),
    )
    fig.update_yaxes(title_text = ATTRIBUTE_COL)

    return dcc.Graph(id="swarm", figure=fig)


@app.callback(
    Output("seq-logos-container", "children"),
    [Input("filtered-data", "children"), Input("hit-data", "children"),],
)
def plot_r_script(filtered_data, hit_data):
    """call R script and render plot"""

    if hit_data is None:
        print("not updating logos")
        raise dash.exceptions.PreventUpdate

    df_filtered = pd.read_json(filtered_data, orient="split")
    df_filtered = df_filtered.where(pd.notnull(df_filtered), None)
    df_hit = pd.read_json(hit_data, orient="split")
    df_hit = df_hit.where(pd.notnull(df_hit), None)

    if len(df_filtered) == len(df_hit):
        raise dash.exceptions.PreventUpdate

    df_not_hit = df_filtered[df_filtered['variant'].isin(set(df_filtered['variant']) - set(df_hit['variant']))]

    # make one big string w/ all sequences, comma separated
    hit_seqs = df_hit["sequence"].tolist()
    not_hit_seqs = df_not_hit["sequence"].tolist()

    hit_seqs = ",".join(hit_seqs)
    not_hit_seqs = ",".join(not_hit_seqs)

    ###### DiffLogo ######

    # Not hit set vs. hit set 
    with tempfile.NamedTemporaryFile() as temp_file:
        # calls R wrapper for DiffLogo that accepts 3 command line arguments, saves png in temp file
        subprocess.call(["./DiffLogo.R", temp_file.name, not_hit_seqs, hit_seqs])

        # convert image to string and pack in URI
        data = base64.b64encode(open(temp_file.name, "rb").read())
        data_url_DiffLogo_not_hit = "data:image/png;base64,{}".format(quote(data))
 
    return html.Div([
            html.H5("Motiv differences between non-hit (top) and hit set (bottom) (DiffLogo)"),
            html.Div(
                html.Img(src=data_url_DiffLogo_not_hit),
                style={"overflow": "auto", "paddingRight": "20px", "padding-left": "20px"}
            ),
    ])


if __name__ == "__main__":
    app.run_server(host='127.0.0.1', port=8050, debug=True)
