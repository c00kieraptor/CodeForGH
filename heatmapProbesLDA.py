#regions are words, patients are documents. Finding probes of interest in topics == finding regions of contributing to that topic.


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy
import statistics

####Read tables####

regScrPrtopic = pd.read_csv("RegScrPrtopic_BRCA-US_FULL.csv", delimiter=",")
regScrPrTopicDrpd = regScrPrtopic.drop(["seqnames", "start", "end", "width", "nCounts", "nCells"], axis=1)

# regAssigUnormal = pd.read_csv("RegAssigUnormal_BRCA-US_FULL.csv", delimiter=",")
# TregAssigUnormal = regAssigUnormal.T
# TregAssigUnormal = TregAssigUnormal.iloc[1: , :]

#Extract first col of regScrPrtopicDrpd
dfForHM = pd.DataFrame(columns = list(regScrPrTopicDrpd))

for COLUMN in regScrPrTopicDrpd.columns[1:]:
    Median = statistics.median(regScrPrTopicDrpd[COLUMN])
    Mean = statistics.mean(regScrPrTopicDrpd[COLUMN])
    Max = max(regScrPrTopicDrpd[COLUMN])
    Cutoff = 0.8
    regScrSupProc = regScrPrTopicDrpd[regScrPrTopicDrpd[COLUMN] > Cutoff]                           #Keep only rows with value OVER Cutoff
    dfForHM = pd.concat([dfForHM, regScrSupProc]).drop_duplicates().sort_values("Unnamed: 0")       #Create df with columns after cutoffs
    name1 = COLUMN + "_Raw.png"                                                                     #
    g = sns.violinplot(y = regScrPrTopicDrpd[COLUMN])
    fig1 = g.get_figure()
    fig1.savefig(name1)
    plt.close()
    name2 = COLUMN+"_Processed.png"
    h = sns.violinplot(y = regScrSupProc[COLUMN])
    fig2 = h.get_figure()
    fig2.savefig(name2)
    plt.close()

dfForHM.to_csv("ProbeTopicScore.csv", index=True)


"""
regScrPrTopicDrpd[(regScrPrTopicDrpd[COLUMN]> Cutoff)]







    lstReduced = list(filter(lambda a: a < lstMedian, lstColValues))

    ax = sns.heatmap(lstReduced)

    g = sns.violinplot(x=list)


    print(regScrPrTopicDrpd[column].median())
    print()

#print(scipy.median_absolute_deviation(regScrPrTopicDrpd))
list = regScrPrTopicDrpd[1].tolist()

regScrSupMedian = regScrPrTopicDrpd.drop(regScrPrTopicDrpd.loc[regScrPrTopicDrpd[COLUMN] >= Median].index)
############# TEST AREA:::
Median = statistics.median(regScrPrTopicDrpd["Scores_Topic1"])
regScrSupMedian = regScrPrTopicDrpd.drop(regScrPrTopicDrpd.loc[regScrPrTopicDrpd["Scores_Topic1"] < Median].index)
g = sns.violinplot(x=regScrSupMedian["Unnamed: 0"], y = regScrSupMedian["Scores_Topic1"])


#regScrPrTopicDrpd.loc[regScrPrTopicDrpd[COLUMN] < Cutoff]
"""