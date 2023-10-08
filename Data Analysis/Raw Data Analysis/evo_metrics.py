import pandas as pd
import progressbar
from pathlib import Path
import final_metrics
import utilities
import csv
import evostatistics

def make_evo_data_dict(filepath, depth):
    data_dict = {}
    metric = utilities.get_metric_from_hjson(filepath)
    key = metric
    if str(filepath).split("\\")[-1] == "EVO_Performance.csv":
        value = make_evo_perf_dict(filepath)
    else:
        value = get_value_from_evo_csv(filepath)
        data_dict.update({"Evo-Epoch":int(str(filepath).split("_")[-1].strip(".csv"))})
    data_dict.update({key:value})
    data_dict = utilities.make_general_dict_from_settings(data_dict, filepath)
    return data_dict

def make_evo_perf_dict(filepath):
    perf_dict ={}
    df = pd.read_csv(filepath)
    for i in range(df.shape[0]):
        perf_dict.update({int(df.at[i,"File:"].split("_")[-1].strip(".pop")):float(df.at[i,"Metric_Score"])})
    return perf_dict

def get_value_from_evo_csv(filepath):
    with open(filepath) as csv_file:
        reader = csv.reader(csv_file)
        try:
            for line in reader:
                value = line[0]
            value = float(value)
        except:
            print("Something went wrong while reading the following file {0}".formate(filepath))
            print("The found value is: {}".formate(value))
            print("Please check this file carefully!")
    return value

def make_df_from_evos(Evopathlist, test=True, ignore_test_range=True):
    EvoDf = pd.DataFrame()
    mainpath = Path.cwd()

    # Stuff for progressbar
    max = len(Evopathlist)
    if test == True:
        test_str = "test"
    else:
        test_str = "train"
    print("Making dataframe for Evo-Steps {0}-performance".format(test_str))
    bar = progressbar.ProgressBar(max_value=max)
    counter = 0

    for path in Evopathlist:
        counter += 1
        depth = final_metrics.get_pathdepth(mainpath, path)
        R2Dict = make_evo_data_dict(path, depth)
        if test == True:
            R2Dict.update({"Set": "Test"})
        else:
            R2Dict.update({"Set": "Train"})
        EvoDf = EvoDf.append(R2Dict, ignore_index=True)
        true_metric = utilities.get_metric_from_hjson(path)
        bar.update(counter)
    if test == True:
        EvoDf = deconvolute_Evo_dict(EvoDf, true_metric)
    EvoDf = utilities.remove_equal_cols(EvoDf, true_metric, ignore_test_range=ignore_test_range)
    print("")
    print("The metric was: {0}".format(true_metric))
    print("")
    print("")
    return EvoDf

def deconvolute_Evo_dict(EvoDf, metric):
    dictlist = []
    for i in range(EvoDf.shape[0]):
        stemdict = EvoDf.iloc[i].to_dict()
        stemdict.pop(str(metric))
        metric_dict = EvoDf.at[i, str(metric)]
        for key in metric_dict:
            evodict = {}
            evodict.update(stemdict)
            evodict.update({"Evo-Epoch": int(key)})
            evodict.update({str(metric): float(metric_dict[key])})
            dictlist.append(evodict)
    newDf = pd.DataFrame(dictlist)
    return newDf

def get_evo_train_values(ignore_test_range=True):
    evoTrainlist = sorted(Path.cwd().glob('**/*EVOFP/EVOstep*.csv'))
    evo_train_df = make_df_from_evos(evoTrainlist, test=False, ignore_test_range=ignore_test_range)
    return evo_train_df

def get_evo_test_values(ignore_test_range=True):
    evoTestlist = sorted(Path.cwd().glob('**/*EVOFP/*EVO_Performance.csv'))
    evo_test_df = pd.DataFrame()
    evo_test_df = make_df_from_evos(evoTestlist, test=True, ignore_test_range=ignore_test_range)
    return evo_test_df

def get_final_evo_values(metric, test=True, train=True, sort=True, save=True, sort_evolast=True, ignore_test_range=True):
    train_df = pd.DataFrame()
    test_df = pd.DataFrame()
    if train:
        train_df = get_evo_train_values(ignore_test_range=ignore_test_range)
    if test:
        test_df = get_evo_test_values(ignore_test_range=ignore_test_range)
    final_df = train_df.append(test_df)
    if sort:
        final_df = evostatistics.order_col_dataframe(final_df, metric, Evolast=sort_evolast)
        final_df = evostatistics.sort_dataframe(final_df, metric)
        final_df = final_df.reset_index(drop=True)
    final_df.fillna(0, inplace=True)
    if save:
        output = "evo_" + str(metric) + ".xlsx"
        final_df.to_excel(output)
    return final_df
