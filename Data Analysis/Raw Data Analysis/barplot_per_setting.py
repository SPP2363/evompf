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

def check_if_col_same(settings_dict, i, df, metrics, wanted_metric, stat_metrics, stats):
    settings_dict = deepcopy(settings_dict)
    del_words = []
    del_words.extend(metrics)
    del_words.append("file_path")
    if stats:
        del_words.append(wanted_metric)
        del_words.extend([*stat_metrics])
    del_words = set(del_words)
    for del_word in del_words:
        settings_dict = del_from_settings(settings_dict, del_word)
    for key in settings_dict:
        a=df.at[i,key]
        b=settings_dict[key]
        if df.at[i,key] != settings_dict[key]:
            return False
    return True

def make_dict_from_row(df, i):
    settings_dict = {}
    for col in df.columns:
        settings_dict.update({col:df.at[i, col]})
    return settings_dict

def contract_dataframe(df, metrics, wanted_metric, stat_metrics, stats):
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
        elif not check_if_col_same(settings_dict, i+1, df, metrics, wanted_metric, stat_metrics, stats):
                settings_dict.update(metric_dict)
                conc_list.append(settings_dict)
                settings_dict = make_dict_from_row(df, i+1)
                for metric in metrics:
                    metric_dict.update({str(metric): []})
    conv_df = pd.DataFrame(conc_list)
    return conv_df



def name_setter(row, metric, wanted_metric, stats):
    name = ""
    headers = list(row.index)
    if stats:
        undesired = ["file_path", "testset_end", "explained_variance_score", "max_error", "mean_squared_error", "mean_absolute_error", "testset_start", "r2_score", "average", "median", "highest", "lowest", "std_deviation", "varriance", "confidence_upp", "confidence_low"]
    else:
        undesired = ["file_path", "testset_end", "explained_variance_score", "max_error", "mean_squared_error", "mean_absolute_error", "r2_score", "average", "median", "highest", "lowest", "std_deviation", "varriance", "confidence_upp", "confidence_low"]
    undesired.append(metric)
    for header in headers:
        if header not in undesired:
            name += str(row[header]) + "_"
    name += metric + "_"
    name += wanted_metric
    name += "_bar.png"
    return name

def build_plot_path(row, metric, wanted_metric, stats):
    name = name_setter(row, metric, wanted_metric, stats)
    if stats:
        save_path = Path.cwd().joinpath("plots", "bar_per_set_stats", wanted_metric)
    else:
        save_path = Path.cwd().joinpath("plots", "bar_per_set", wanted_metric)
    if not save_path.exists():
        Path.mkdir(save_path, parents=True)
    file = save_path.joinpath(name)
    return file

def sorter(other_settings, input_df, stats):
    sort_by = other_settings
    sort_by.append("Set")
    if "testset_start" in list(input_df.columns.values):
        sort_by.append("testset_start")
    if "sheet_name" in list(input_df.columns.values):
        sort_by.append("sheet_name")
    if len(sort_by) > 0:
        input_df_sorted = input_df.sort_values(sort_by, ignore_index=True)
    else:
        input_df_sorted = input_df.copy(deep=True)
    return input_df_sorted

def format_label(wanted_metric):
    if wanted_metric == "r2_score":
        return "R$^2$"
    if wanted_metric == "mean_absolute_error":
        return "MAE"
    if wanted_metric == "mean_squared_error":
        return "MSE"
    if wanted_metric == "max_error":
        return "Max Error"
    if wanted_metric == "explained_variance_score":
        return "/u03C3$^2_{exp}$"
    if wanted_metric == "feature_type":
        return "feature type"
    if wanted_metric == "fp_size":
        return "fingerprint size"
    if wanted_metric == "newgen_rate":
        return "mutation rate"

def bar_plotter(row, metric, wanted_metric, stats, width_bar=0.2, color="royalblue", set_ax_manually=False, yaxis_lowlim=0, yaxis_uplim=1, fontsize=11, **kwargs):
    file = build_plot_path(row, metric, wanted_metric, stats)
    x = range(len(row[metric]))
    xaxis_lowlim = min(x) - 0.5
    xaxis_uplim = max(x) + 0.5
    if stats:
        height = row["average"]
        lolims = row["confidence_low"]
        uplims = row["confidence_upp"]
        if Plot_Objects.bar_per_set.error:
            error = [(uplims[i] - lolims[i]) / 2 for i in range(len(lolims))]
        else:
            error = None
    else:
        height = row[wanted_metric]
        error = None
    fig, ax = plt.subplots()
    ax.bar(x, height, yerr=error, align="center", width=width_bar, color=color, error_kw={"capsize":4.0, "capthick":1.5, "ecolor":"black", "elinewidth":1.5})
    ax.set_xticks(x)
    ax.set_xticklabels(row[metric], fontsize=fontsize)
    ax.set_ylabel(format_label(wanted_metric), fontsize=fontsize)
    ax.set_xlabel(format_label(metric), fontsize=fontsize)
    ax.set_xlim([xaxis_lowlim, xaxis_uplim])
    if set_ax_manually:
        ax.set_ylim([yaxis_lowlim, yaxis_uplim])
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_axisbelow(True)
    #plt.show()
    plt.savefig(file, format="png", bbox_inches="tight")
    plt.close()


def bar_plot_manager(input_df, wanted_metric, stats, stat_metrics=None, **kwargs):
    settings = Plot_Objects.bar_per_set.settings
    for setting in settings:
        if setting in list(input_df.columns.values):
            other_settings = [item for item in Plot_Objects.bar_per_set.settings.copy() if item in list(input_df.columns.values)]
            other_settings.remove(setting)
            input_df_sorted = sorter(other_settings, input_df, stats)
            if stats:
                metrics = [setting, "average", "confidence_upp", "confidence_low"]
            else:
                metrics = [setting, wanted_metric]
            df = contract_dataframe(input_df_sorted, metrics, wanted_metric, stat_metrics, stats)
            for i, row in df.iterrows():
                bar_plotter(row, setting, wanted_metric, stats, **kwargs)