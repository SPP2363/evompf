import sys
from pathlib import Path
import pandas as pd
from rdkit import Chem
from statistics import variance, mean
from progressbar import progressbar
from math import sqrt
from rdkit.Chem.Draw import rdMolDraw2D, MolToImage
from IPython.display import SVG, display
from rdkit.Chem import Draw



#SMARTS einlesen
smarts_dataset = pd.read_excel("SMARTS_Evo.xlsx")
smarts_dataset["SMARTS"] = smarts_dataset["Query"].apply(lambda x: Chem.MolFromSmarts(x))

#DD dataset to molecules
dd_dataset = pd.read_excel("Dreher_and_Doyle_input_data.xlsx")
dd_dataset["Ligand_mol"] = dd_dataset["Ligand"].apply(lambda x: Chem.MolFromSmiles(x))
dd_dataset["Additive_mol"] = dd_dataset["Additive"].apply(lambda x: Chem.MolFromSmiles(x))
dd_dataset["Base_mol"] = dd_dataset["Base"].apply(lambda x: Chem.MolFromSmiles(x))
dd_dataset["Aryl halide_mol"] = dd_dataset["Aryl halide"].apply(lambda x: Chem.MolFromSmiles(x))

def draw_molecules_and_smarts(mol, smarts):
    Draw.MolToFile(mol,"Test.png")



#Visualize Mols/SMARTS
for i in progressbar(range(smarts_dataset.shape[0])):
    smarts = smarts_dataset.at[i,"SMARTS"]
    for j in range(dd_dataset.shape[0]):
        ligand = dd_dataset.at[j, "Ligand_mol"]
        additive = dd_dataset.at[j, "Additive_mol"]
        base = dd_dataset.at[j, "Base_mol"]
        halide = dd_dataset.at[j, "Aryl halide_mol"]

        ligand_checker = ligand.HasSubstructMatch(smarts)
        if ligand_checker:
            smiles = Chem.MolToSmiles(dd_dataset.at[j, "Ligand_mol"])
            name = (str(smarts_dataset.at[i, "Query"]) + "_" + str(smiles) + ".png").replace(":","_")
            path = Path.cwd().joinpath(name)
            Draw.MolToFile(ligand, path)

        additive_checker = additive.HasSubstructMatch(smarts)
        if additive_checker:
            smiles = Chem.MolToSmiles(dd_dataset.at[j, "Additive_mol"])
            name = (str(smarts_dataset.at[i, "Query"]) + "_" + str(smiles) + ".png").replace(":","_")
            path = Path.cwd().joinpath(name)
            Draw.MolToFile(additive, path)

        base_checker = base.HasSubstructMatch(smarts)
        if base_checker:
            smiles = Chem.MolToSmiles(dd_dataset.at[j, "Base_mol"])
            name = (str(smarts_dataset.at[i, "Query"]) + "_" + str(smiles) + ".png").replace(":","_")
            path = Path.cwd().joinpath(name)
            Draw.MolToFile(base, path)

        halide_checker = halide.HasSubstructMatch(smarts)
        if halide_checker:
            smiles = Chem.MolToSmiles(dd_dataset.at[j, "Aryl halide_mol"])
            name = (str(smarts_dataset.at[i, "Query"]) + "_" + str(smiles) + ".png").replace(":","_")
            path = Path.cwd().joinpath(name)
            Draw.MolToFile(halide, path)



#Match SMARTS on molecules & calc avg. yield
pos_means = []
pos_variances = []
neg_means = []
neg_variances = []
nr_pos = []
nr_neg = []
for i in progressbar(range(smarts_dataset.shape[0])):
    pos_yields = []
    neg_yields = []
    smarts = smarts_dataset.at[i,"SMARTS"]
    for j in range(dd_dataset.shape[0]):
        ligand = dd_dataset.at[j, "Ligand_mol"]
        additive = dd_dataset.at[j, "Additive_mol"]
        base = dd_dataset.at[j, "Base_mol"]
        halide = dd_dataset.at[j, "Aryl halide_mol"]
        ligand_checker = ligand.HasSubstructMatch(smarts)
        additive_checker = additive.HasSubstructMatch(smarts)
        base_checker = base.HasSubstructMatch(smarts)
        halide_checker = halide.HasSubstructMatch(smarts)
        if ligand_checker or additive_checker or base_checker or halide_checker:
            rxn_yield = float(dd_dataset.at[j, "Output"])
            pos_yields.append(rxn_yield)
        else:
            rxn_yield = float(dd_dataset.at[j, "Output"])
            neg_yields.append(rxn_yield)
    pos_variance_value =sqrt(variance(pos_yields))
    pos_variances.append(pos_variance_value)
    pos_mean_value = mean(pos_yields)
    pos_means.append(pos_mean_value)
    neg_variance_value =sqrt(variance(neg_yields))
    neg_variances.append(neg_variance_value)
    neg_mean_value = mean(neg_yields)
    neg_means.append(neg_mean_value)
    nr_pos.append(len(pos_yields))
    nr_neg.append(len(neg_yields))
smarts_dataset["average_yield_pos"] = pos_means
smarts_dataset["average_yield_neg"] = neg_means
smarts_dataset["std_yield_pos"] = pos_variances
smarts_dataset["std_yield_neg"] = neg_variances
smarts_dataset["nr_datapoints_match"] = nr_pos
smarts_dataset["nr_datapoints_unmatch"] = nr_neg
smarts_dataset.to_excel("Yield_averadge_and_STD.xlsx")
print("")