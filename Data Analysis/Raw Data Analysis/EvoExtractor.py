import os
import re
import pandas as pd
from progressbar import progressbar
from pathlib import Path

#in subordner ausführen
MainP = os.getcwd()
MainP = os.path.join(MainP, "Data")
MainP = os.path.join(MainP,"DD_Mutation_0.25_Lencheck")
os.chdir(MainP)

class Scores():
    R2test: float
    R2train: float
    Setting: str
    Split: str
    Experiment: str
    Path: Path

    def __init__(self):
        self.R2test = 0.0
        self.R2train = 0.0
        self.Setting = ""
        self.Split = ""
        self.Experiment = ""
        self.Path = os.getcwd()

    def __str__(self):
        return "R2test: " + str(self.R2test) + "\n" + "R2train: " + str(self.R2train) + "\n" "Path: " + str(
            self.Path) + "\n" "Settings: " + str(self.Setting) + "\n" "Split: " + str(
            self.Split) + "\n" "Experiment: " + str(self.Experiment)


class EvoScore():
    Score: Scores  # R2train will be R2 of periode
    step: int

    def __init__(self):
        self.Score = Scores()
        self.step = 0

    def __str__(self):
        return "Periode: " + str(self.step) + "\n" + str(self.Score)


def getR2(pathoffile, EvoOnOff=False):
    if EvoOnOff is False:
        df = pd.read_csv(pathoffile)
        return df.at[0, "r2_score"]
    else:
        df = pd.read_csv(pathoffile)

        return df.at[df.index[-1],"FP_0"]


def compareScores(Score1, Score2):
    if Score1.Split == Score2.Split:
        if Score1.Experiment == Score1.Experiment:
            return True
        else:
            print("Experiment was not equal")
            return False
    else:
        print("Splitting was not equal")
        return False


def setScore(PathOfTestCsv, PathOfTrainCsv=None, EvoOnOff=False):
    score = Scores()

    score.R2test = getR2(PathOfTestCsv, EvoOnOff)
    if PathOfTrainCsv is not None:
        score.R2train = getR2(PathOfTrainCsv, EvoOnOff)
    if EvoOnOff is False:
        if PathOfTestCsv.parents[1] == PathOfTrainCsv.parents[1]:
            score.Path = PathOfTestCsv
            score.Setting = PathOfTestCsv.parents[1].name
            score.Split = PathOfTestCsv.parents[2].name
            score.Experiment = PathOfTestCsv.parents[3].name
        else:
            print("WARNING: Paths of test and train unequal - check if file is missing")
            print("Return empty score class object")
            score = Scores()
            return score
    else:
        score.Path = PathOfTestCsv
        score.Setting = PathOfTestCsv.parents[2].name
        score.Split = PathOfTestCsv.parents[3].name
        score.Experiment = PathOfTestCsv.parents[4].name

    return score


def setEvoScore(PathOfEvoCSV):
    EvoS = EvoScore()
    EvoS.Score = setScore(PathOfEvoCSV, None, True)
    EvoS.step = (int(re.findall(r'\d+',PathOfEvoCSV.name)[0]))
    return EvoS


def setListofEvoScores(PathOfEvoPerformance): #Takes "EVO-Performance.csv and makes List of EvoScore-Objects out of it
    df = pd.read_csv(PathOfEvoPerformance, index_col=False)
    EvoS = EvoScore()
    ListofEvoScores = []
    for i in range(df.shape[0]):
        EvoS = EvoScore()
        EvoS.Score.Path = Path(df.at[i,"File:"])
        EvoS.Score.Setting = EvoS.Score.Path.parents[2].name
        EvoS.Score.Split = EvoS.Score.Path.parents[3].name
        EvoS.Score.Experiment = EvoS.Score.Path.parents[4].name
        EvoS.Score.R2test = df.at[i,"Metric_Score"]
        EvoS.step = (int(re.findall(r'\d+',EvoS.Score.Path.name)[0]))
        ListofEvoScores.append(EvoS)
    return ListofEvoScores


def scoretodict(Scores):
    scoreDict = {"Setting": Scores.Setting, "Split": Scores.Split, "Experiment": Scores.Experiment,
                 "R2-Train": Scores.R2train, "R2-Test": Scores.R2test}
    return scoreDict


def evotodict(Evo):
    evoDict = {"Step": Evo.step, "Setting": Evo.Score.Setting, "Split": Evo.Score.Split, "Experiment": Evo.Score.Experiment,
                 "R2-Test": Evo.Score.R2test}
    return evoDict


# Inizitizing
R2 = False
Evo = True

settingsStrToInt = True
evoS = EvoScore()

if R2 is True:
    Trainlist = sorted(Path.cwd().glob('**/*Trainprediction.csv'))
    Testlist = sorted(Path.cwd().glob('**/*Testprediction.csv'))
    print("Making AllScores.csv-file with all R2-Scores.")
    print("")
    R2Datalist: list = []
    for i, path in enumerate(Trainlist):  # Enumerate gibt zwei Werte zurück 0: index 1:element der list
        R2Datalist.append(setScore(Testlist[i], path))
    AuxDF = pd.DataFrame(columns=["Experiment", "Split", "Setting", "R2-Train", "R2-Test"])
    for score in R2Datalist:
        AuxDF = AuxDF.append(scoretodict(score), ignore_index=True)
    with open("AllScores.csv", "w") as file:
        file.write("Path , " + str(Trainlist[0].parents[3]) + "\n")
        file.write("\n")
        for i in AuxDF.Experiment.unique():
            file.write(str(i) + "\n")
            file.write("\n")
            ExDF = AuxDF.loc[AuxDF["Experiment"] == i]
            for j in AuxDF.Split.unique():
                file.write(str((j) + "\n"))
                file.write("Setting , R2-Train , R2-Test" + "\n")
                ExSplitDF = ExDF.loc[AuxDF["Split"] == j]
                settingsStrToInt = True
                ExSplitDF = ExSplitDF.reset_index(drop=True)
                for l in range(ExSplitDF.shape[0]):
                    try:
                        ExSplitDF.at[l, "Setting"] = int(ExSplitDF.at[l, "Setting"])
                    except:
                        settingsStrToInt = False
                if settingsStrToInt is True:
                    ExSplitDF = ExSplitDF.sort_values("Setting")
                ExSplitDF = ExSplitDF.reset_index(drop=True)
                for k in range(ExSplitDF.shape[0]):
                    file.write(str(
                        str((ExSplitDF.at[k, "Setting"])) + " , " + str(ExSplitDF.at[k, "R2-Train"]) + " , " + str(
                            ExSplitDF.at[k, "R2-Test"]) + "\n"))
                file.write("\n")
            file.write("\n")
else:
    print("Analysation of R2's is turned off")
    print("")

if Evo is True:

    # Part for getting EVO_Performance.csv to nice List of EvoScore Elements
    R2TestEvoList: list = []
    R2TestPathList = []
    R2TestPathList = sorted(Path.cwd().glob('**/EVO_Performance.csv'))
    EvoAuxDFtest = pd.DataFrame(columns=["Experiment", "Split", "Setting", "Step", "R2-Test"])
    for i, element in enumerate(R2TestPathList):
        R2TestEvoList = R2TestEvoList + setListofEvoScores(element)


    print("Making Dicts from your R2-Test-Scores")
    for i, element in enumerate(R2TestEvoList):
        EvoAuxDFtest = EvoAuxDFtest.append(evotodict(element), ignore_index=True)
    print("Done")
    print("")

    #Part for getting Files with Train-R2 in Dataframe
    print("Making Evo-Scores.csv-file with all Evo-Scores.")
    print("")
    EvoDatalist: list = []
    Evolist = sorted(Path.cwd().glob('**/EVOstep*.csv'))
    print("Listing all Files:")
    for i in progressbar(range(len(Evolist))):  # Enumerate gibt zwei Werte zurück 0: index 1:element der list
        EvoDatalist.append(setEvoScore(Evolist[i]))
    print(len(EvoDatalist), "Files were found.")
    print("")
    EvoAuxDF = pd.DataFrame(columns=["Experiment", "Split", "Setting", "Step", "R2-Test"])
    print("Generating dataframe-object:")
    for evo in progressbar(EvoDatalist):
        EvoAuxDF = EvoAuxDF.append(evotodict(evo), ignore_index=True)
    print("Generated")
    print("")
    print("Writing to Evo-Scores.csv-file:")
    EvoAuxDF = EvoAuxDF.rename(columns ={"R2-Test":"R2-Train"}, inplace=False)
    EvoAuxDF = EvoAuxDF.reindex(columns = EvoAuxDF.columns.tolist()+ ["R2-Test"])
    print("")
    print("Merging dataframes")
    for i in progressbar(range(EvoAuxDF.shape[0])):
        for j in range(EvoAuxDFtest.shape[0]):
            if str(EvoAuxDF.at[i,"Experiment"]) == str(EvoAuxDFtest.at[j,"Experiment"]) and str(EvoAuxDF.at[i,"Split"]) == str(EvoAuxDFtest.at[j,"Split"]) and str(EvoAuxDF.at[i,"Setting"]) == str(EvoAuxDFtest.at[j,"Setting"]) and str(EvoAuxDF.at[i,"Step"]) == str(EvoAuxDFtest.at[j,"Step"]):
                EvoAuxDF.at[i, "R2-Test"] = EvoAuxDFtest.at[j, "R2-Test"]

    print("")
    print("Done")
    print("")

    with open("EvoScores.csv", "w") as file:
        file.write("Path , " + str(Evolist[0].parents[3]) + "\n")
        file.write("\n")
        for i in progressbar(EvoAuxDF.Experiment.unique()):
            file.write(str(i) + "\n")
            file.write("\n")
            EvoExDF = EvoAuxDF.loc[EvoAuxDF["Experiment"] == i]
            for j in EvoExDF.Split.unique():
                file.write(str((j) + "\n"))
                EvoExSplitDF = EvoExDF.loc[EvoExDF["Split"] == j]
                settingList = EvoExSplitDF.Setting.unique()
                settingsStrToInt = True
                for k in range(len(settingList)):
                    try:
                        settingList[k] = int(settingList[k])
                    except:
                        settingsStrToInt = False
                if settingsStrToInt is True:
                    settingList.sort()
                printDF = pd.DataFrame()
                maxlength = 0
                sortdfList = []
                dfDict = {}
                sortDict = {}
                #Erzeugen einer Liste die alle zu schreibenden Elemente enthält - Column name ist dabei gewählte Einstellung
                for l, element in enumerate(settingList):
                    preprintDF = EvoExSplitDF.loc[EvoExSplitDF["Setting"] == str(element)].sort_values("Step")
                    preprintDF = preprintDF[["Setting","Step","R2-Train","R2-Test"]].reset_index(drop=True) #Wähle bestimmte Column aus
                    preprintDF = preprintDF.rename(columns={"Step":str(element)}, inplace=False) #Nennt Column Step um in Setting
                    preprintDF = preprintDF[[preprintDF.at[0, "Setting"],"R2-Train","R2-Test"]]
                    if preprintDF.shape[0] > maxlength:
                        maxlength = preprintDF.shape[0]
                    dfDict.update({str(element):preprintDF})
                sortdfList = sorted(dfDict, reverse=True)

                #List häufig Str -> Alle, wenn möglich zu ints machen
                for l, stri in enumerate(sortdfList):
                    try:
                        sortdfList[l] = int(stri)
                    except:
                        str = str
                sortdfList = sorted(sortdfList)
                for l, inte in enumerate(sortdfList):
                    sortdfList[l] = str(inte)

                #Schleife darüber wie häufig geschrieben werden muss (Maximale EvoSteps)
                for key in sortdfList:
                    file.write(str(dfDict.get(key).columns[0]) + " , R2-Train , R2-Test , ")
                file.write("\n")
                for m in range(maxlength):
                    RowList = []
                    # Schleife um sich aus jedem DF der dFListe die Elemente zu ziehen
                    for key in sortdfList:
                        columns = dfDict.get(key).columns
                        try:
                            RowList.append(dfDict.get(key).at[m, columns[0]])
                            RowList.append(dfDict.get(key).at[m, columns[1]])
                            RowList.append(dfDict.get(key).at[m, columns[2]])
                        except:
                            RowList.append("")
                            RowList.append("")
                            RowList.append("")
                    for n in range(len(RowList)):
                        file.write(str(RowList[n]) + " , ")
                    file.write("\n")
                file.write("\n")
                file.write("\n")
            file.write("\n")
        file.write("\n")
        print("All data is now easily available for you.")
               # file.write("Setting , R2-Train , R2-Test" + "\n")
               # ExSplitDF = ExDF.loc[AuxDF["Split"] == j]
               # settingsStrToInt = True
               # ExSplitDF = ExSplitDF.reset_index(drop=True)
else:
    print("Analysation of Evo-Scores is turned off.")
    print("")
