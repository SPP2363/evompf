from pathlib import Path
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from statistics import variance, mean
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

df_path = Path.cwd().joinpath("Pattern_investigation.xlsx")
if not df_path.is_file():
    """Count Patterns"""

    #Read EVO-

    EVO_files = list(Path.cwd().glob("*/*.csv"))
    cols = [f"FP_{i}" for i in range(80)]
    evopattern_string = ""

    for EVO_file in EVO_files:
        dataset = pd.read_csv(EVO_file, sep=",",usecols=cols)
        evopattern_string += dataset.to_string()

    atomic_numbers = [i+1 for i in range(118)]
    nr_of_pos_atom_patterns_per_atom = [] #print_to_dataset
    nr_of_neg_atom_patterns_per_atom = [] #print_to_dataset

    #Count per Atom
    for atomic_number in atomic_numbers:
        nr = atomic_number
        reg_positive = f"[^!]#{nr}[^0-9]"
        reg_negative = f"!#{nr}[^0-9]"
        list_of_pos_atom_patterns = re.findall(reg_positive,evopattern_string)
        list_of_neg_atom_patterns = re.findall(reg_negative,evopattern_string)
        nr_of_pos_atom_patterns = len(list_of_pos_atom_patterns)
        nr_of_neg_atom_patterns = len(list_of_neg_atom_patterns)
        nr_of_pos_atom_patterns_per_atom.append(nr_of_pos_atom_patterns)
        nr_of_neg_atom_patterns_per_atom.append(nr_of_neg_atom_patterns)

    """Count number of samples with atom and tox statistics"""
    tox_path = Path.cwd().joinpath("ToxTraining.xlsx")
    tox_dataset = pd.read_excel(tox_path)
    SMARTS_str_list = [f"[#{i+1}]" for i in range(118)]
    SMARTS_list = [Chem.MolFromSmarts(smarts_string) for smarts_string in SMARTS_str_list]
    nr_of_molecules_with_atom = [0]*118 #print_to_dataset

    LD_50_values_per_atom = []
    LogLD_50_values_per_atom = []
    for i in range(118):
        LD_50_values_per_atom.append([])
        LogLD_50_values_per_atom.append([])

    for i in range(tox_dataset.shape[0]):
        smiles = tox_dataset.at[i,"SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logLD50 = tox_dataset.at[i,"LD50"]
        LD50 = (10**logLD50)*mw
        for i, SMARTS in enumerate(SMARTS_list):
            hit_list = mol.GetSubstructMatches(SMARTS)
            if hit_list:
                nr_of_molecules_with_atom[i] = nr_of_molecules_with_atom[i]+1
                LD_50_values_per_atom[i].append(LD50)
                LogLD_50_values_per_atom[i].append(logLD50)

    """Statistical data (mean/variance)"""
    LD_50_mean_list = [0]*118
    logLD_50_mean_list = [0]*118
    LD_50_variance_list = [0]*118
    logLD_50_variance_list = [0]*118

    for i in range(118):
        LD_50_values = LD_50_values_per_atom[i]
        LogLD_50_values = LogLD_50_values_per_atom[i]
        if len(LogLD_50_values) > 1:

            LD_50_mean = mean(LD_50_values)
            LogLD_50_mean = mean(LogLD_50_values)
            LD_50_variance = variance(LD_50_values)
            LogLD_50_variance = variance(LogLD_50_values)

            LD_50_mean_list[i] = LD_50_mean
            logLD_50_mean_list[i] = LogLD_50_mean
            LD_50_variance_list[i] = LD_50_variance
            logLD_50_variance_list[i] = LogLD_50_variance


    cols_new_df = ["Element","Nr positive Evo Patterns","Nr negative Evo Patterns","Number of SMILES","LD50 mean", "logLD50 mean","LD50 variance", "logLD50 variance"]
    df_content = list(zip(SMARTS_str_list,nr_of_pos_atom_patterns_per_atom,nr_of_neg_atom_patterns_per_atom,nr_of_molecules_with_atom,LD_50_mean_list,logLD_50_mean_list,LD_50_variance_list,logLD_50_variance_list))
    final_dataframeset = pd.DataFrame(df_content, columns=cols_new_df)
    final_dataframeset.to_excel(df_path)

final_dataframeset = pd.read_excel(df_path)

def log_and_normalize_colum(df, old_colum, new_colum):
    col_list = []
    for i in range(df.shape[0]):
        value = df.at[i,old_colum]
        if value < 1:
            value = 1
        logvalue = math.log10(value)
        col_list.append(logvalue)
    df[new_colum] = col_list
    min_v = df[new_colum].min()
    max_v = df[new_colum].max()
    df[new_colum] = (df[new_colum]-min_v)/(max_v-min_v)
    return df

#Drop Element with occurence in Dataset below n_smiles
n_smiles = 3
droplist = []
n_int = np.nan
for i in range(final_dataframeset.shape[0]):
    if final_dataframeset.at[i,"Number of SMILES"] < n_smiles:
        final_dataframeset.at[i, 'Nr positive Evo Patterns'] = 0
        final_dataframeset.at[i, 'logLD50 mean'] = float("Nan")
        final_dataframeset.at[i, 'logLD50 variance'] = float("Nan")
        final_dataframeset.at[i, "Number of SMILES"] = 0

#Normalize Values for colour plots -> SMILES and nr of patterns are in log
final_dataframeset["logLD50 mean norm"] = (((final_dataframeset["logLD50 mean"]-final_dataframeset["logLD50 mean"].min())/(final_dataframeset["logLD50 mean"].max()-final_dataframeset["logLD50 mean"].min())-1)*-1)
final_dataframeset["logLD50 variance norm"] = (final_dataframeset["logLD50 variance"]-final_dataframeset["logLD50 variance"].min())/(final_dataframeset["logLD50 variance"].max()-final_dataframeset["logLD50 variance"].min())
final_dataframeset = log_and_normalize_colum(final_dataframeset, "Number of SMILES", "Number of SMILES norm")
final_dataframeset = log_and_normalize_colum(final_dataframeset, "Nr positive Evo Patterns", "Nr positive Evo Patterns norm")

#Crazy PSE Plot
elements = {1:"H",2:"He",3:"Li",4:"Be",5:"B",6:"C",7:"N",8:"O",9:"F",10:"Ne",11:"Na",12:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar",19:"K",20:"Ca",21:"Sc",22:"Ti",23:"V",24:"Cr",
            25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",31:"Ga",32:"Ge",33:"As",34:"Se",35:"Br",36:"Kr",37:"Rb",38:"Sr",39:"Y",40:"Zr",41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",
            47:"Ag",48:"Cd",49:"In",50:"Sn",51:"Sb",52:"Te",53:"I",54:"Xe",55:"Cs",56:"Ba",57:"La",72:"Hf",73:"Ta",74:"W",75:"Re",76:"Os",77:"Ir",78:"Pt",79:"Au",80:"Hg",81:"Tl",
            82:"Pd",83:"Bi",84:"Po",85:"At",86:"Rn"}

PSE_mapping = {1:[1,1],2:[1,18]
            ,3:[2,1],4:[2,2],5:[2,13],6:[2,14],7:[2,15],8:[2,16],9:[2,17],10:[2,18],
            11:[3,1],12:[3,2],13:[3,13],14:[3,14],15:[3,15],16:[3,16],17:[3,17],18:[3,18],
            19:[4,1],20:[4,2],21:[4,3],22:[4,4],23:[4,5],24:[4,6],25:[4,7],26:[4,8],27:[4,9],28:[4,10],29:[4,11],30:[4,12],31:[4,13],32:[4,14],33:[4,15],34:[4,16],35:[4,17],36:[4,18],
            37:[5,1],38:[5,2],39:[5,3],40:[5,4],41:[5,5],42:[5,6],43:[5,7],44:[5,8],45:[5,9],46:[5,10],47:[5,11],48:[5,12],49:[5,13],50:[5,14],51:[5,15],52:[5,16],53:[5,17],54:[5,18],
            55:[6,1],56:[6,2],57:[6,3],72:[6,4],73:[6,5],74:[6,6],75:[6,7],76:[6,8],77:[6,9],78:[6,10],79:[6,11],80:[6,12],81:[6,13],82:[6,14],83:[6,15],84:[6,16],85:[6,17],86:[6,18],
            }

nRows = 6
nCols = 18
PSE_dataframe = pd.DataFrame(index=range(nRows),columns=range(nCols))

#Collect Data per Element
for element_number in PSE_mapping:
    element_symbol = elements.get(element_number)
    coordinates = PSE_mapping.get(element_number)
    if coordinates:
        row = coordinates[0]-1 #For dicts i counted like a human -> -1 here corrects it
        col = coordinates[1]-1
        nr_of_pos_atom_patterns_element = final_dataframeset.at[element_number-1,"Nr positive Evo Patterns"]
        nr_of_smiles_element = final_dataframeset.at[element_number - 1, "Number of SMILES"]
        logLD50_mean_element = final_dataframeset.at[element_number - 1, "logLD50 mean"]
        logLD50_var_element = final_dataframeset.at[element_number - 1, "logLD50 variance"]
        nr_of_pos_atom_patterns_element_norm = final_dataframeset.at[element_number - 1, "Nr positive Evo Patterns norm"]
        nr_of_smiles_element_norm = final_dataframeset.at[element_number - 1, "Number of SMILES norm"]
        logLD50_mean_element_norm = final_dataframeset.at[element_number - 1, "logLD50 mean norm"]
        logLD50_var_element_norm = final_dataframeset.at[element_number - 1, "logLD50 variance norm"]
        data_dict = {"symbol":element_symbol,"number":element_number,"nr_of_pos_atom_patterns_element":nr_of_pos_atom_patterns_element,"nr_of_smiles_element":nr_of_smiles_element,
                    "logLD50_mean_element":logLD50_mean_element,"logLD50_var_element":logLD50_var_element,"nr_of_pos_atom_patterns_element norm":nr_of_pos_atom_patterns_element_norm,
                     "nr_of_smiles_element norm":nr_of_smiles_element_norm, "logLD50_mean_element norm":logLD50_mean_element_norm,"logLD50_var_element norm":logLD50_var_element_norm}
        PSE_dataframe.at[row,col] = data_dict

#Expand DF
expended_df = pd.DataFrame(index=range(nRows*2),columns=range(nCols*2),dtype=float)
non_data_dict = {"symbol": float("Nan"), "number": float("Nan"),
             "nr_of_pos_atom_patterns_element": float("Nan"),
             "nr_of_smiles_element": float("Nan"),
             "logLD50_mean_element": float("Nan"), "logLD50_var_element": float("Nan"),
             "nr_of_pos_atom_patterns_element norm": float("Nan"),
             "nr_of_smiles_element norm": float("Nan"),
             "logLD50_mean_element norm": float("Nan"),
             "logLD50_var_element norm": float("Nan")}

for i in range(nRows):
    for j in range(nCols):
        data_dict = PSE_dataframe.at[i,j]
        if not isinstance(data_dict, dict):
            data_dict = non_data_dict
        expended_df.at[i * 2, j * 2] = float(data_dict["nr_of_pos_atom_patterns_element norm"])
        expended_df.at[i * 2 + 1, j * 2] = float(data_dict["nr_of_smiles_element norm"])
        expended_df.at[i * 2, j * 2 + 1] = float(data_dict["logLD50_mean_element norm"])
        expended_df.at[i * 2 + 1, j * 2 + 1] = float(data_dict["logLD50_var_element norm"])

fig, ax = plt.subplots()
fig.set_size_inches(10, 3.3)

#color_bar
N = 256
name = 'epoch_cmap'
col1 = (128, 0, 0)
col2 = (0, 42, 88)
col3= (19, 79, 13)
step1=0.25
step2=0.25
step_fade=0.05
col_lst = [col1, col2, col3]
num_col = len([x for x in col_lst if x is not None])
color_inter = np.ones((int(step1 * N), 4))
color_inter[:, 0] = np.linspace(col1[0] / 256, col2[0] / 256, int(step1 * N))
color_inter[:, 1] = np.linspace(col1[1] / 256, col2[1] / 256, int(step1 * N))
color_inter[:, 2] = np.linspace(col1[2] / 256, col2[2] / 256, int(step1 * N))
color_inter2 = np.ones((int(step2 * N), 4))
color_inter2[:, 0] = np.linspace(col2[0] / 256, col3[0] / 256, int(step2 * N))
color_inter2[:, 1] = np.linspace(col2[1] / 256, col3[1] / 256, int(step2 * N))
color_inter2[:, 2] = np.linspace(col2[2] / 256, col3[2] / 256, int(step2 * N))
color_inter = np.vstack((color_inter, color_inter2))
cm.register_cmap(name=name, cmap=ListedColormap(np.vstack(color_inter)))

#pltsettings
im = ax.pcolor(expended_df, cmap = 'epoch_cmap',zorder=0)
#plt.rcParams["figure.figsize"] = (50,5)

#axis and frame
ax.invert_yaxis()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

#lines between elements
x, y, w, h = 0.0, 1, 18, 1
for key in PSE_mapping:
    x = (PSE_mapping[key][1]-1)
    y = (PSE_mapping[key][0]-1)
    ax.add_patch(Rectangle((x*2, y*2), 2, 2, fill=False, edgecolor='lightgrey', lw=1.5, clip_on=False, zorder=5))

#add PSE-Symbols
for i in range(nRows):
    for j in range(nCols):
        data_dict = PSE_dataframe.at[i,j]
        if isinstance(data_dict, dict):
            top_left = data_dict["symbol"]
            text = ax.text(j * 2+1,i * 2+1,top_left,ha="center", va="center", color="w",fontsize=16, fontweight="bold",zorder=6,fontfamily='Century Gothic')

#hide elements with too little datapoints
ax.add_patch(Rectangle((0, 0), 2, 2, fill=True, facecolor='grey', edgecolor='lightgrey', lw=1.5, clip_on=False,zorder=3))
for key in PSE_mapping:
    y = (PSE_mapping[key][1]-1)
    x = (PSE_mapping[key][0]-1)
    element_info =PSE_dataframe.at[x,y]
    if isinstance(element_info, dict):
        nr_of_smiles_element = element_info['nr_of_smiles_element']
        if nr_of_smiles_element < 3:
            ax.add_patch(Rectangle((y*2, x*2), 2, 2, fill=True, facecolor='grey', edgecolor='lightgrey', lw=1.5, clip_on=False,zorder=3))

#add PSE-Numbers
#add numbers
for i in range(nRows):
    for j in range(nCols):
        data_dict = PSE_dataframe.at[i,j]
        if isinstance(data_dict, dict):
            top_left = data_dict["nr_of_pos_atom_patterns_element"]
            position_y_top_left = (j*2)+0.5
            position_x_top_left = (i*2)+0.25
            text = ax.text(position_y_top_left, position_x_top_left, top_left, ha="center", va="center", color="w", fontsize=6,zorder=2,fontfamily='Century Gothic')

            top_right = round(data_dict["logLD50_mean_element"],2)
            position_y_top_right = (j * 2) + 1.5
            position_x_top_right = (i * 2) + 0.25
            ax.text(position_y_top_right, position_x_top_right, top_right, ha="center", va="center", color="w", fontsize=6,zorder=2,fontfamily='Century Gothic')

            bottom_left = round(data_dict["nr_of_smiles_element"],0)
            position_y_bottom_left = (j * 2) + 0.5
            position_x_bottom_left = (i * 2) + 1.75
            ax.text(position_y_bottom_left,position_x_bottom_left , bottom_left, ha="center", va="center", color="w", fontsize=6,zorder=2,fontfamily='Century Gothic')

            bottom_right = round(data_dict["logLD50_var_element"],2)
            position_y_bottom_right = (j * 2) + 1.5
            position_x_bottom_right = (i * 2) + 1.75
            ax.text(position_y_bottom_right,position_x_bottom_right, bottom_right, ha="center", va="center", color="w", fontsize=6,zorder=2,fontfamily='Century Gothic')

fig.tight_layout()
my_dpi = 300
plt.savefig("PSE_plot_long.png", format="png", bbox_inches="tight", transparent=True, dpi=my_dpi)
#plt.show()
