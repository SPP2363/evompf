from pathlib import Path
import pandas as pd
import hjson
from copy import deepcopy
from collections import OrderedDict

def getParentName(pathoffile:Path, level:int):
    name = str(pathoffile.parents[level]).split("\\")[-1]
    return name

def converte_ints_in_list(list:list):
    for i in range(len(list)):
        try:
            list[i] = int(list[i])
        except:
            list[i] = list[i]
        return list

def check_if_keyword_in_list(list:list, keyword:str):
    checker = False
    for element in list:
        if keyword in element:
            checker = True
        elif keyword == element:
            checker = True
    return checker

def flatter_settingsdict(settings_dict, datadict = {}):
    auxdict = deepcopy(settings_dict)
    for key in auxdict:
        if isinstance(auxdict[key], OrderedDict) and bool(auxdict[key]):
            datadict.update(flatter_settingsdict(auxdict[key], datadict))
        elif bool(auxdict[key]):
            newkey = key
            counter = 0
            while newkey in datadict:
                counter += 1
                newkey = str(newkey) + str("_") + str(counter)
            datadict.update({newkey:auxdict[key]})
    return datadict

def make_general_dict_from_settings(data_dict, filepath):
    if "EVO" in filepath.name and "prediction" not in filepath.name:
        jsonpath = filepath.parents[2].joinpath("settings.hjson")
    else:
        jsonpath = filepath.parents[1].joinpath("settings.hjson")
    with open(jsonpath, "r") as jason_file:
        settings_dict = hjson.load(jason_file)
    data_dict.update(flatter_settingsdict(settings_dict, {}))
    return data_dict

def make_general_data_dict(data_dict, filepath, depth):
    for i in range(depth):
        parentname = getParentName(filepath,i+1)
        if "_" in parentname:
            parent_list = parentname.split("_")
            parent_list = converte_ints_in_list(parent_list)

            #check for other files which have underscores, ints in first item or keywords (Iter/iter/Depth/depth/cat/Cat/RF/rf/
            keywordlist_aglo = ["Iter", "iter", "Depth", "depth", "cat", "Cat", "RF", "rf", "CB", "cb", "10K", "50K","default"]
            keywordlist_test = ["test","Test"]
            checker_algo = False
            checker_test = False
            for keyword in keywordlist_test:
                checker_test = check_if_keyword_in_list(parent_list, keyword)
                if checker_test == True:
                    break
            for keyword in keywordlist_aglo:
                checker_algo = check_if_keyword_in_list(parent_list, keyword)
                if checker_algo == True:
                    break
            if checker_test:
                value = parentname
                key = "splitting"
            elif isinstance(parent_list[0], int) or checker_algo:
                value = parentname
                key = "algo-settings"
            else:
                value = parent_list[-1]
                del parent_list[-1]
                key = '_'.join(parent_list)

        else:
            test = True
            try:
                int(parentname)
            except:
                test = False
            if test:
                key = "Number"
                value = parentname
            else:
                key = str(i)
                value = parentname
        data_dict.update({key:value})
    return data_dict

def remove_equal_cols(df, metric, ignore_test_range=True):
    collist = df.columns
    for column in collist:
        try:
            if column == "Set":
                pass
            elif column == metric:
                pass
            elif column == "thread_count" or column == "thread_count_1":
                df.drop([column], axis=1, inplace=True)
            elif df[column].nunique() == 1:
                df.drop([column], axis=1, inplace=True)
        except:
            df.drop([column], axis=1, inplace=True)
    if "sheet_name_1" in df.columns:
        df.drop(["sheet_name_1"], axis=1, inplace=True)
    collist = df.columns
    if "sheet_name" in collist:
        if ignore_test_range:
            if "testset_start" in collist:
                df.drop(["testset_start"], axis=1, inplace=True)
            if "testset_end" in collist:
                df.drop(["testset_end"], axis=1, inplace=True)
    return df

def check_incomplete(Evo=False):
    hjson_list = sorted(Path.cwd().glob('**/settings.hjson'))
    wronglist = []
    print("Checking for incomplete calculations")
    for path in hjson_list:
        perf_check = True
        train_check = True
        test_check = True
        perf_path_a = path.parents[0].joinpath("Output").joinpath("CV_EVOFP").joinpath("EVO_Performance.csv")
        perf_path_b = path.parents[0].joinpath("Output").joinpath("CV0_EVOFP").joinpath("EVO_Performance.csv")
        if len(sorted(path.parents[0].glob('**/*Trainprediction.csv'))) != 1:
            train_check = False
        if len(sorted(path.parents[0].glob('**/*Testprediction.csv'))) != 1:
            test_check = False
        if not perf_path_b.is_file() and not perf_path_a.is_file():
            perf_check = False
        if Evo == False:
            perf_check = True
        if perf_check == False or train_check == False or test_check == False:
            wronglist.append({"Missing_file": str(path), "Performance": perf_check, "Test_csv": test_check, "Train_csv": train_check})
    if wronglist:
        print("WARNING!!!! - Incomplete calulation found. Please consider the Missingfiles.xlsx")
        wrong_df = pd.DataFrame(wronglist)
        wrong_df.to_excel("Missingfiles.xlsx")
    else:
        wrong_df = pd.DataFrame()
        print("Cool! All experiments look good!")
    print("")
    print("")
    return wrong_df

def get_metric_from_hjson(filepath):
    jsonpath = filepath.parents[2].joinpath("settings.hjson")
    with open(jsonpath, "r") as jason_file:
        settings_dict = hjson.load(jason_file)
    if settings_dict["encoder_settings"]["encoder_settings"]["fitfunc_name"] == "diversity":
        metric = "diversity"
    else:
        metric = settings_dict["encoder_settings"]["encoder_settings"]["fitfunc_param"]["metricscore"]
    return metric

def get_evo_metric():
    evoTestlist = sorted(Path.cwd().glob('**/*EVOFP/*EVO_Performance.csv'))
    metric_set= set()
    for filepath in evoTestlist:
        metric = get_metric_from_hjson(filepath)
        metric_set.add(metric)
    if len(metric_set) < 1:
            print("WARNING: Evometrics seem to be not the same. Found metrics: {}".format(metric_set))
            print("The following metric will be used: {}".format(metric))
    return metric
