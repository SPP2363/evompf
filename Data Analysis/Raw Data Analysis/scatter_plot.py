import os
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from plots_class import *


def testprediction_reader():
    testlist = sorted(Path.cwd().glob('**/*Testprediction.csv'))
    test = True
    return testlist, test

def trainprediction_reader():
    trainlist = sorted(Path.cwd().glob('**/*Trainprediction.csv'))
    test = False
    return trainlist, test

def name_setter(csv_file, test):
    path = csv_file.parents[1]
    depth = Plot_Objects.scatter.depth_name
    name = str(path.relative_to(csv_file.parents[depth+1]))
    print()
    name = name.replace("\\", "-")
    name += "-scatter"
    if test:
        name += "_test.png"
    if not test:
        name += "_train.png"
    return name

def create_xy_table(csv_file):
    df = pd.read_csv(csv_file)
    xy_table = df.filter(["True_values", "Predicted_values"])
    return xy_table

def build_plot_path(csv_file, test):
    name = name_setter(csv_file, test)
    if test:
        set = "test"
    else:
        set = "train"
    save_path = Path.cwd().joinpath("plots", "scatter", set)
    if not save_path.exists():
        Path.mkdir(save_path, parents=True)

    file = save_path.joinpath(name)
    return file

def create_line(xy_table, ax):
    x_line = [xy_table["True_values"].min(), xy_table["True_values"].max()]
    y_line = [xy_table["True_values"].min(), xy_table["True_values"].max()]
    ax.plot(x_line, y_line, '--', color="grey")

def scatter_plotter(csv_file, test, xaxis_lowlim = 0, xaxis_uplim = 10, yaxis_lowlim = 0, yaxis_uplim = 10, fontsize=11, set_ax_manually=False, **kwargs):
    xy_table = create_xy_table(csv_file)
    file = build_plot_path(csv_file, test)
    fig, ax = plt.subplots()
    create_line(xy_table, ax)
    ax.scatter(xy_table["True_values"], xy_table["Predicted_values"], marker=".", s = 40, c = "royalblue", alpha=.5, zorder=2.5)
    ax.set_xlabel(Plot_Objects.scatter.xlabel + "/" + Plot_Objects.scatter.unit, fontsize=fontsize)
    ax.set_ylabel(Plot_Objects.scatter.ylabel + "/" + Plot_Objects.scatter.unit, fontsize=fontsize)
    if set_ax_manually:
        ax.set_xlim([xaxis_lowlim, xaxis_uplim])
        ax.set_ylim([yaxis_lowlim, yaxis_uplim])
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    #plt.grid(b=True, which="both")
    #ax.set_axisbelow(True)
    #plt.show()
    plt.savefig(file, format="png", bbox_inches="tight")
    plt.close()

def scatter_plot_manager(**kwargs):
    if Plot_Objects.scatter.test:
        testlist, test = testprediction_reader()
        for csv_file in testlist:
            scatter_plotter(csv_file, test, **kwargs)
    if Plot_Objects.scatter.train:
        trainlist, test = trainprediction_reader()
        for csv_file in trainlist:
            scatter_plotter(csv_file, test, **kwargs)




#Testarea
if __name__ == "__main__":
    main_path = Path("C:\\Users\\felix\\Documents\\AK\\Research\\evo_fp\\Test_for_plotter\\feature_typ_match")
    os.chdir(main_path)
    scatter_plot_manager(Plot_Objects)
    print()