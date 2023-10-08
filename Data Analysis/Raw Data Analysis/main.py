import os
import hjson
import pandas as pd
from pathlib import Path
import evostatistics
import final_metrics
import evo_metrics
import utilities
import epoch_stats_plot
import epoch_plot
import scatter_plot
import barplot_per_setting
from plots_class import *

    #todo plots
    #todo checker what should be written to file (modular usage)
    #checker call every module individually!


def main():

    jsonpath = Path.cwd().joinpath("config.hjson")
    with open(jsonpath, "r") as jason_file:
        settings= hjson.load(jason_file)

    main_path = settings["path"]
    save = settings["test"]
    test = settings["test"]
    train = settings["train"]
    evo = settings["evo"]
    statistics = settings["stats"]
    sort = settings["sort"]
    plots = settings["plots"]
    check_missing= settings["check_missing"]
    check_similarity= settings["check_similarty"]

    sort_evolast = settings["sort_evolast"]
    wanted_metrics = settings["metric"]
    stat_metrics = settings["stat_metrics"]
    ignore_test_range = settings["ignore_test_range"]

    os.chdir(main_path)

    print("working in path {}".format(str(main_path)))
    print("")
    for wanted_metric in wanted_metrics:
        print("The program is working with {} as metric.".format(wanted_metric))
        print("")
        '''This part of the program will extract all metric values from the corresponding .csv files'''
        if train or test:
            final_metric_df = final_metrics.get_final_values(metric= wanted_metric, test=test, train=train, sort=sort, save=save, ignore_test_range=ignore_test_range)

        '''This part of the program will extract all evo-values from the corresponding .csv files'''
        if evo:
            evo_metric = utilities.get_evo_metric()
            evo_metric_df = evo_metrics.get_final_evo_values(metric = evo_metric, test=test, train=train, sort=sort, save=save, sort_evolast=sort_evolast, ignore_test_range= ignore_test_range)
        else:
            evo_metric_df = pd.DataFrame()

        '''This part of the program will check if there are settings-files without .csvs'''
        if check_missing:
            wronglist = utilities.check_incomplete(Evo=evo)

        '''programm part for statistics, if multiple experiments were made'''
        if statistics:
            stats_dict, evo_stats_dict = evostatistics.all_stat_metrics(final_metric_df, wanted_metric, stat_metrics=stat_metrics,  evo=evo, evo_df=evo_metric_df, check_settings_similarity= check_similarity, ignore_test_range=ignore_test_range)

        if plots:
            if Plot_Objects.scatter.draw == True:
                kwargs = Plot_Objects.scatter.optional_settings
                scatter_plot.scatter_plot_manager(**kwargs)
            if Plot_Objects.epoch_stats.draw == True and statistics:
                if wanted_metric == "r2_score":
                    kwargs = Plot_Objects.epoch_stats.optional_settings
                    epoch_stats_plot.epoch_stats_plot_manager(evo_stats_dict, wanted_metric, stat_metrics, **kwargs)
            if Plot_Objects.epoch.draw == True:
                if evo and wanted_metric == "r2_score":
                    kwargs = Plot_Objects.epoch.optional_settings
                    epoch_plot.epoch_plot_manager(evo_metric_df, wanted_metric, ignore_test_range, **kwargs)
            if Plot_Objects.bar_per_set.draw == True:
                kwargs = Plot_Objects.bar_per_set.optional_settings
                if statistics:
                    stats_df = stats_dict["stats"]
                    barplot_per_setting.bar_plot_manager(stats_df, wanted_metric, statistics, stat_metrics=stat_metrics, **kwargs)
                else:
                    barplot_per_setting.bar_plot_manager(final_metric_df, wanted_metric, statistics, **kwargs)




    #todo übergeben manuelle values einstellungen für alte jobs
if __name__ == "__main__":
    main()