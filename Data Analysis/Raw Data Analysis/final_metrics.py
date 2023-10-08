import pandas as pd
import progressbar
from pathlib import Path
import evostatistics
import utilities

def check_metric(pathoffile:Path, metric):
    columns = pd.read_csv(pathoffile).columns
    if metric in columns:
        return True
    else:
        return False

def check_metric_pathlist(pathlist, metric):
    stop = False
    if not pathlist:
        stop = True
        print("Error: There are no searchable experiments!")
    for path in pathlist:
        if not check_metric(path, metric):
            print("Error: Your Metric \"{}\" does not exist in the output-file".format(metric))
            stop = True
            break
    return stop

def get_pathdepth(mainpath, filepath):
    #Depth marks the number of extra columns you need in the dataframe
    depth = 0
    while True:
        if str(filepath.parents[depth]) == str(mainpath):
            return depth-1
        depth += 1

def get_value(pathoffile:Path, metric:str):
    df = pd.read_csv(pathoffile)
    metric_value = float(df.at[0, metric])
    return metric_value

def make_data_dict(filepath, depth, metric):
    data_dict = {}
    key = metric
    value = get_value(filepath, metric)
    data_dict.update({key:value})
    utilities.make_general_dict_from_settings(data_dict, filepath)
    #utilities.make_general_data_dict(data_dict, filepath, depth)
    return data_dict

def make_df_from_final_csv(Pathlist, metric, test=True, ignore_test_range=True):
    #Stuff for progressbar
    mainpath = Path.cwd()
    R2Df = pd.DataFrame()
    max = len(Pathlist)
    if test == True:
        test_str = "test"
    else:
        test_str = "train"
    print("Making dataframe the following metric: {1} on the {0}-set".format(test_str, str(metric)))
    #bar = progressbar.ProgressBar(max_value=max)
    counter = 0

    for path in Pathlist:
        counter += 1
        depth = get_pathdepth(mainpath, path)
        R2Dict = make_data_dict(path, depth, metric)
        if test == True:
            R2Dict.update({"Set": "Test"})
        else:
            R2Dict.update({"Set": "Train"})
        R2Dict.update({"file_path":str(path)})
        R2Df = R2Df.append(R2Dict, ignore_index=True)
        #bar.update(counter)

    #Aussortieren columns mit nur gleichen Einträgen
    R2Df = utilities.remove_equal_cols(R2Df, metric, ignore_test_range=ignore_test_range)
    print("")
    print("")
    return R2Df

def get_train_values(metric, ignore_test_range=True):
    trainlist = sorted(Path.cwd().glob('**/*Trainprediction.csv'))
    stop_train = check_metric_pathlist(trainlist, metric)
    if not stop_train:
        R2Df = make_df_from_final_csv(trainlist, metric, test=False, ignore_test_range=ignore_test_range)
    else:
        R2Df = pd.DataFrame()
    return R2Df

def get_test_values(metric, ignore_test_range=True):
    testlist = sorted(Path.cwd().glob('**/*Testprediction.csv'))
    stop_train = check_metric_pathlist(testlist, metric)
    if not stop_train:
        R2Df = make_df_from_final_csv(testlist, metric, test=True, ignore_test_range=ignore_test_range)
    else:
        R2Df = pd.DataFrame()
    return R2Df

def get_final_values(metric, test=True, train=True, sort= True, save=True, ignore_test_range=True):
    train_df = pd.DataFrame()
    test_df = pd.DataFrame()
    if train:
        train_df = get_train_values(metric, ignore_test_range=ignore_test_range)
    if test:
        test_df = get_test_values(metric, ignore_test_range=ignore_test_range)
    final_df = train_df.append(test_df)
    if sort:
        final_df = evostatistics.order_col_dataframe(final_df, metric)
        final_df = evostatistics.sort_dataframe(final_df, metric)
        final_df = final_df.reset_index(drop=True)
    final_df.fillna(0, inplace=True)
    if save:
        output = "final_" + str(metric) + ".xlsx"
        final_df.to_excel(output)
    return final_df


