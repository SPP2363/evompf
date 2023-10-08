import os
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from plots_class import *
from copy import deepcopy

def del_from_settings(settings_dict, del_word):
    if del_word in settings_dict:
        del settings_dict[del_word]
    return settings_dict

def check_if_col_same(settings_dict, i, df, metrics, stat_metrics, wanted_metric, ignore_test_range=True):
    settings_dict = deepcopy(settings_dict)
    del_words = []
    del_words.extend(metrics)
    del_words.append("file_path")
    del_words.append(wanted_metric)
    del_words.extend([*stat_metrics])
    if ignore_test_range:
        del_words.append("testset_end")
        del_words.append("testset_start")
    for del_word in del_words:
        settings_dict = del_from_settings(settings_dict, del_word)
    for key in settings_dict:
        if df.at[i,key] != settings_dict[key]:
            return False
    return True

def make_dict_from_row(df, i):
    settings_dict = {}
    for col in df.columns:
        settings_dict.update({col:df.at[i, col]})
    return settings_dict

def contract_dataframe(df, metrics, stat_metrics, wanted_metric, ignore_test_range= True):
    conc_list = []
    metric_dict = {}
    for metric in metrics:
        metric_dict.update({str(metric):[]})
    settings_dict = make_dict_from_row(df,0)
    df_len = df.shape[0]
    for i in range(df_len):
        for metric in metrics:
            metric_dict[metric].append(df.at[i, metric])
        if i == df_len-1:
            settings_dict.update(metric_dict)
            conc_list.append(settings_dict)
        elif not check_if_col_same(settings_dict, i+1, df, metrics, stat_metrics, wanted_metric, ignore_test_range=ignore_test_range):
                settings_dict.update(metric_dict)
                conc_list.append(settings_dict)
                settings_dict = make_dict_from_row(df, i+1)
                for metric in metrics:
                    metric_dict.update({str(metric): []})
    conv_df = pd.DataFrame(conc_list)
    return conv_df



def name_setter(row):
    name = ""
    headers = list(row.index)
    undesired = ["testset_start", "testset_end", "Evo-Epoch", "r2_score", "average", "median", "highest", "lowest", "std_deviation", "varriance", "confidence_upp", "confidence_low"]
    for header in headers:
        if header not in undesired:
            name += str(row[header]) + "_"
    name += "epoch_stats.png"
    return name

def build_plot_path(row):
    name = name_setter(row)
    save_path = Path.cwd().joinpath("plots", "epoch_stats")
    if not save_path.exists():
        Path.mkdir(save_path, parents=True)
    file = save_path.joinpath(name)
    return file

def epoch_plotter(row, bars=True, set_ax_manually=False, xaxis_lowlim = 0, xaxis_uplim = 52, yaxis_lowlim = 0, yaxis_uplim = 1, fontsize=11, **kwargs):
    file = build_plot_path(row)
    x = row["Evo-Epoch"]
    y = row["average"]
    lolims = row["confidence_low"]
    uplims = row["confidence_upp"]
    error = [(uplims[i] - lolims[i])/2 for i in range(len(lolims))]
    if bars:
        fig, ax = plt.subplots()
        ax.errorbar(x, y, yerr=error, marker='o', ms=5, mec="mediumblue", mfc="lightblue", linestyle='none', capsize=3.0, capthick=0.6, ecolor="black", elinewidth=0.7)
    else:
        fig, ax = plt.subplots()
        plt.scatter(x, y, alpha=.8, c=error, cmap='seismic')
        plt.colorbar(label="Confidence 90%")
    ax.set_xlabel(Plot_Objects.epoch_stats.xlabel, fontsize=fontsize)
    ax.set_ylabel(Plot_Objects.epoch_stats.ylabel, fontsize=fontsize)
    if set_ax_manually:
        ax.set_xlim([xaxis_lowlim, xaxis_uplim])
        ax.set_ylim([yaxis_lowlim, yaxis_uplim])
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.grid(visible=True, which="both")
    ax.set_axisbelow(True)
    #plt.show()
    plt.savefig(file, format="png", bbox_inches="tight")
    plt.close()


def epoch_stats_plot_manager(evo_stats_dict, wanted_metric, stat_metrics, **kwargs):
    metrics = Plot_Objects.epoch_stats.metrics
    df = contract_dataframe(evo_stats_dict["evo_stats"], metrics, stat_metrics, wanted_metric, ignore_test_range = True)
    for i, row in df.iterrows():
        epoch_plotter(row, **kwargs)

#Testarea
'''excel_file = Path("C:\\Users\\felix\\Documents\\AK\\Research\\evo_fp\\Test_for_plotter\\evo_r2_score_stats.xlsx")
main_path = Path("C:\\Users\\felix\\Documents\\AK\\Research\\evo_fp\\Test_for_plotter")
os.chdir(main_path)

df = pd.read_excel(excel_file)
epoch_plot_manager(excel_file,Plot_Objects)
print()'''