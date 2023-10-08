import pandas as pd
import statistics
import math
import scipy.stats as st
import utilities
from copy import deepcopy

def get_ordered_col_list(df, metric, Evolast=True):
    collist = list(df.columns)
    if "Evo-Epoch" in collist:
        collist.remove("Evo-Epoch")
        if Evolast == False:
            collist.insert(0, "Evo-Epoch")
        else:
            collist.append("Evo-Epoch")
    if "Number" in collist:
        collist.remove("Number")
        collist.insert(0, "Number")
    if str(metric) in collist:
        collist.remove(str(metric))
        collist.append(str(metric))
    return collist

def order_col_dataframe(df, metric, Evolast=True):
    new_df = df.copy()
    #evolast=False will sort dataframe in this order that the epochs will be sorted with top one priority (good for epoch-plots)
    collist = get_ordered_col_list(new_df, metric, Evolast)
    new_df = new_df[collist]
    return new_df

def key_to_collist_end(collist, key):
    if key in collist:
        collist.remove(key)
        collist.append(key)
    return collist

def sort_dataframe(df, metric):
    collist = list(df.columns)
    collist = key_to_collist_end(collist, "sheet_name")
    collist = key_to_collist_end(collist, "Set")
    collist.reverse()
    collist = key_to_collist_end(collist,"testset_start")
    collist = key_to_collist_end(collist,"testset_end")
    try:
        collist.remove(str(metric))
        collist.append(str(metric))
    except Exception as e:
        print(e)
    df = df[collist]
    for col in collist:
        if isinstance(df.at[0,col], list):
            collist.remove(col)
    if "file_path" in collist:
        collist.remove("file_path")
    df = df.sort_values(by=collist, ascending=False)
    df.reset_index(drop=True, inplace=True)
    return df

def get_n_samples(df:pd.DataFrame):
    #get number of samples per experiment
    collist = list(df.columns)
    if "Number" in collist:
        df = df.astype({"Number": "int32"})
        n = len(set(df["Number"].tolist()))
        return n
    else:
        return 1

def check_similarity(df, conc_df, metric, n, n_exp):
    if "Number" in list(conc_df.columns):
        test_df = df.drop(columns=["Number"])
    else:
        test_df = df.copy(deep=True)
    test_conc_df = conc_df.drop(columns=[metric])
    test_df.drop(columns=[metric], inplace=True)
    checker = True
    for i in range(n_exp):
        for j in range(n):
            for column in list(test_conc_df.columns):
                if str(test_df.at[(i*(n)+j),column]) != str(test_conc_df.at[i,column]):
                    checker = False
    if checker:
        return True
    else:
        return False

def del_from_settings(settings_dict, del_word):
    if del_word in settings_dict:
        del settings_dict[del_word]
    return settings_dict

def check_if_col_same(settings_dict, i, df, metric, ignore_test_range=True):
    settings_dict = deepcopy(settings_dict)
    del_words = [metric]
    del_words.append("file_path")
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
    settingsdict = {}
    for col in df.columns:
        settingsdict.update({col:df.at[i, col]})
    return settingsdict

def contract_dataframe(df, metric:str, ignore_test_range= True):
    df = order_col_dataframe(df, metric, Evolast=True)
    df = sort_dataframe(df, metric)
    conc_list = []
    metric_list = []
    settings_dict = make_dict_from_row(df,0)
    df_len = df.shape[0]
    for i in range(df_len):
        metric_list.append(df.at[i, metric])
        if i == df_len-1:
            settings_dict.update({metric: metric_list})
            conc_list.append(settings_dict)
        elif not check_if_col_same(settings_dict, i+1, df, metric, ignore_test_range=ignore_test_range):
                settings_dict.update({metric: metric_list})
                conc_list.append(settings_dict)
                settings_dict = make_dict_from_row(df, i+1)
                metric_list = []
    conv_df = pd.DataFrame(conc_list)
    conv_df = sort_dataframe(conv_df, metric)
    return conv_df

def get_confidence(values, percentage=0.95):
    mean = statistics.mean(values)
    std = statistics.stdev(values)
    n = len(values)
    z = st.norm.ppf(((1-percentage)/2)+percentage)
    confidence_up = mean + z * std/math.sqrt(n)
    confidence_low = mean - z * std/math.sqrt(n)
    return  confidence_up, confidence_low

def metric_from_list(values:list, metrics:dict):
    values = [float(element) for element in values]
    metric_dict = {}
    if len(values) == 1:
        return metric_dict
    for metric in metrics:
        if metric == "average":
            metric_dict.update({"average":statistics.mean(values)})
        elif metric == "median":
            metric_dict.update({"median":statistics.median(values)})
        elif metric == "highest":
            metric_dict.update({"highest":max(values)})
        elif metric == "lowest":
            metric_dict.update({"lowest":min(values)})
        elif metric == "std_deviation":
            metric_dict.update({"std_deviation":statistics.stdev(values)})
        elif metric == "varriance":
            metric_dict.update({"varriance":statistics.variance(values)})
        elif metric == "confidence":
            conf_percentage = float(metrics[metric]["percentage"])
            up, low = get_confidence(values, conf_percentage)
            metric_dict.update({"confidence_upp":up})
            metric_dict.update({"confidence_low": low})
        else:
            try:
                methode = getattr(statistics, str(metric))
                metric_value = methode(values)
                metric_dict.update({"metric":metric_value})
            except:
                pass
    return metric_dict

def calc_stat_metrics(df, stat_metrics:dict, metric:str, evo_step= False, check_settings_similarity= True, ignore_test_range=True):
    conc_df = contract_dataframe(df, metric,ignore_test_range= ignore_test_range)
    metric_list = []
    for i in range(conc_df.shape[0]):
        values = conc_df.at[i, metric]
        metric_dict = metric_from_list(values, stat_metrics)
        metric_list.append(metric_dict)
    aux_df = pd.DataFrame(metric_list)
    for column in list(aux_df.columns):
        conc_df[column] = aux_df[column]
    return conc_df

def all_stat_metrics(final_df, metric, stat_metrics:dict,  evo=False, evo_df=pd.DataFrame(),check_settings_similarity= True, ignore_test_range=True):
    print("Okay it seems that you have repeated some experiments!")
    print("Making files with statistical metrics")
    print("")
    print("The following metrics will be included:")
    print(stat_metrics.keys())
    print("")
    if evo:
        print("You choosed to analyse final- and evo-step values")
    else:
        print("You choosed to analyse final values only")
    metric_df = calc_stat_metrics(final_df, stat_metrics, metric, evo_step=False, check_settings_similarity=check_settings_similarity,  ignore_test_range= ignore_test_range)
    stats_output_name = "final_" + str(metric) + "_stats.xlsx"
    metric_df.reset_index(drop=True, inplace=True)
    metric_df.to_excel(stats_output_name)
    if evo:
        evo_metric = utilities.get_evo_metric()
        evo_metric_df = calc_stat_metrics(evo_df, stat_metrics, evo_metric, evo_step=True, check_settings_similarity=check_settings_similarity,  ignore_test_range= ignore_test_range)
        evo_stats_output_name = "evo_" + str(metric) + "_stats.xlsx"
        evo_metric_df.reset_index(drop=True, inplace=True)
        evo_metric_df.to_excel(evo_stats_output_name)
    else:
        evo_metric_df = pd.DataFrame()
    evo_metric_dict = {"final_evo_stats": metric,"evo_stats":evo_metric_df}
    metric_dict = {"final_stats": metric,"stats":metric_df}
    print("")
    print("")
    return metric_dict, evo_metric_dict

def calc_full_predc_metrics():
    pass