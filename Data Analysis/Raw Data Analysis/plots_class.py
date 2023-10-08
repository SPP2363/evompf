from pathlib import Path
import hjson

plot_settings_path = Path.cwd().joinpath("plot_settings.hjson")
with open(plot_settings_path, "r") as plot_settings_file:
    plot_settings= hjson.load(plot_settings_file)

class Plot:
    def __init__(self, dictionary= {}):
        for k, v in dictionary.items():
            setattr(self, k, v)

class Plots:
    plots_dict:dict
    scatter:Plot
    epoch:Plot

    def __init__(self, dictionary):
        self.construct("scatter", dictionary)
        self.construct("epoch", dictionary)
        self.construct("epoch_stats", dictionary)
        self.construct("bar_per_set", dictionary)

    def construct(self, name, dict):
        if dict[name]:
            setattr(self, name, Plot(dict[name]))
        else:
            setattr(self, name, Plot({"draw":False}))

    def append(self, dict):
        key = dict.keys[0]
        self.key = dict[key]


Plot_Objects = Plots(plot_settings)

print()

