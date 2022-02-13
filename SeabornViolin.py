"""
The methylation level of CpGs in the topics, you could plot the complete distribution of values for a given set of patients associated with the corresponding topic and compare to other CpGs in the other topics. For instance, using violin plots with strips
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import lzma
#####


####Read tables#######################
"""
methTable = pd.read_csv("MethylTable_BRCA-US_FULL.csv.xz", delimiter=",", compression="xz")
regAssigUnormal = pd.read_csv("RegAssigUnormal_BRCA-US_FULL.csv", delimiter=",")
regScrPrtopic = pd.read_csv("RegScrPrtopic_BRCA-US_FULL.csv", delimiter=",")
topicAssigToPatient = pd.read_csv("topicAssigToPatient_BRCA-US_FULL.csv", delimiter=",")

"""
methTable = pdf.read_csv(snakemake.input[methTab], delimiter=",", compression="xz")
regScrPrtopic = pd.read_csv(snakemake.input[regScrNorm], delimiter=",")
regAssigUnormal = pd.read_csv(snakemake.input[regScrUnrm], delimiter=",")
topicAssigToPatient = pd.read_csv(snakemake.input[topicAssig], delimiter=",")

####################################


"""
Following to find the donors that have the highest topic assignment:
"""
#Remove unnamed column with topic numbers and set those as rownames
topicAsPatRownames = topicAssigToPatient.set_index("Unnamed: 0")
topicAsPatRownames.index.names = [None]
#topicAsPatRownames["mean"] = topicAsPatRownames.mean(axis=1)

#subract mean of each 
topicPatMean = topicAsPatRownames.sub(topicAsPatRownames.mean(axis=1), axis=0)
	
#topicAbsMean = topicAbsMean.abs()


#Make dictionary-class .
class topDict(dict):
	def __init__(self):
		self = dict()
	def add(self, col, row):
		self[col] = row	





# Make a list of column names
colList = list(topicPatMean)



topicPatDict = topDict()

#Add key (patientID) with n topics (as values), after subtracting mean of that topic - this because it can find the topics contributing most for that patient compared to other topics.
for item in colList:
	topicPatDict.col = item
	patient = topicPatMean[item]
	topicPatDict.row = patient.nlargest(4)
	topicPatDict.add(topicPatDict.col, list(topicPatDict.row.index))





#Find number of topics 
tpcNo = len(topicAssigToPatient)

#Make lists for each topic, named topic_ + topic number 
nr = 0

while nr < tpcNo: 
	nr = nr + 1
	locals()[f"Topic_{str(nr)}"] = []
	

#number of topics and add 1 because this is used for setting columns to be worked on on following for loop. (First column is region: "Unnamed: 0")
tcpNoPlus1 = tpcNo+1

#Make list with Topic_X as name of list (x = topic number) and patients as items where the patients had that topic as highest contributing topic. This will cluster patients that has the same "highest contributing topic" together
for key, value in topicPatDict.items():
    for n in range(1,tcpNoPlus1,1):
        if n in value:
            locals()[f"Topic_{str(n)}"].append(key)



"""
Following to find regions most contributing to a topic
"""
#Set region as rownumber methTable and drop the column. Also make regTable as DF with only region as column so that columns from topics in next while-loop can be added
methTabWithRegAsRowname = methTable
methTabWithRegAsRowname.index = methTable["Unnamed: 0"] 
methTabWithRegAsRowname.index.names = [None]

regTable = methTabWithRegAsRowname["Unnamed: 0"]

methTabWithRegAsRowname = methTabWithRegAsRowname.drop(["Unnamed: 0"], axis=1)


###################################    RegScrPrTopic    #######################################################
#Remove unneeded columns
regScrPrTopicDrpd = regScrPrtopic.drop(["seqnames", "start", "end", "width", "nCounts", "nCells"], axis=1)

#Make dict where keys are columns(topicScores) and values are columns of "Unnamed: 0 and the values, respectively. For the values, those below 0.5, are discarded "
LookupDict = {}
for col in regScrPrTopicDrpd.columns[1:]:
    newDF = regScrPrTopicDrpd[["Unnamed: 0", col]]
    newDF = newDF[newDF[col] > 0.5]
    LookupDict[col] = list(newDF["Unnamed: 0"])
"""
{'Scores_Topic1':                       Unnamed: 0
								205     chrX:134569361-134569362
								227       chrX:89177497-89177498
								314       chrX:75249221-75249222
								387     chrX:111874417-111874418
								389     chrX:145077290-145077291
"""
#####################################################################################
###


#Make dataframe named TopFrame_X with X_patients as columns to ease further downstream analysis (X = topic number). Make new DF of 30 % of columns in each Topic
nr = 0
while nr < tpcNo:
	nr = nr + 1
	locals()[f"Scores_Topic{str(nr)}"] = methTabWithRegAsRowname[methTabWithRegAsRowname.columns & locals()[f"Topic_{str(nr)}"]]
	locals()[f"Scores_Topic{str(nr)}"] = locals()[f"Scores_Topic{str(nr)}"].add_prefix(f"{str(nr)}_") 
	locals()[f"Scores_Topic{str(nr)}"] = pd.concat([methTable["Unnamed: 0"], (locals()[f"Scores_Topic{str(nr)}"].sample(frac=1, axis=1))], axis=1)


#Remove rows from Scores_TopicX that are not present in the LookupDict
for key in LookupDict.keys():
    locals()[f"{str(key)}"] = (locals()[key].loc[LookupDict[key]])






nr = 0
while nr < tpcNo:
    nr = nr + 1
    df = locals()[f"Scores_Topic{str(nr)}"].round(decimals=2)
    for column in df.columns[1:]:
        df1 = pd.concat([df["Unnamed: 0"], df[column]], axis=1)
        df1 = df1.drop_duplicates(subset=[column])
        g = sns.violinplot(y=column, data=df1)
        fig = g.get_figure()
        name = f"{column}.png"
        fig.savefig(name)
        plt.clf()







#less decimals in the table to easily drop duplicates in the next for-loop
deciRegTab = regTable.round(decimals=3)	 
dRTSansUn = deciRegTab.drop(["Unnamed: 0"], axis=1)

#drop rows that have duplicate values in dataframe. Then make new dataframe of each colum with that column and column named "Unnamed: 0"
for column in dRTSansUn.columns:
	locals()[f"DF_{column}"] = pd.concat([deciRegTab["Unnamed: 0"], dRTSansUn[column]], axis=1)
	locals()[f"DF_{column}"] = locals()[f"DF_{column}"].drop_duplicates(subset=[column])
	g = sns.violinplot(y=column, data=locals()[f"DF_{column}"])
	fig = g.get_figure()
	name=column+".png"
	fig.savefig(name)
	plt.clf()

"""	
#drop rows that have duplicate values in dataframe. Then make new dataframe of each colum with that column and column named "Unnamed: 0"
for column in dRTSansUn.columns:
	locals()[f"tab_{column}"] = pd.concat([deciRegTab["Unnamed: 0"], dRTSansUn[column]], axis=1)
	locals()[f"tab_{column}"] = locals()[f"tab_{column}"].drop_duplicates(subset=[column])
	g = sns.violinplot(y=column, data=locals()[f"tab_{column}"])
	fig = g.get_figure()
	name=column+".png"
	fig.savefig(name)
	plt.clf()
"""


	


#g = sns.violinplot(y="1_DO2395", data=tab_1_DO2395)
"""
#Add last column with average value of every row
regTabWAv = regTable
regTabWAv['Average'] = regTabWAv.mean(axis=1)


#Find difference between rows and average
diffFrame = regTabWAv.diff(axis=1, periods=-1)

#Find lowest value in row
minValDiffF = diffFrame.min(axis=1)

#Round decimals and remove dup values for each column
#BUT SHOULD REMOVE SAME FOR ALL WITH DUPE THERE!!



#Violin plot 

df = regTable.melt("Unnamed: 0", var_name="Patient",  value_name="Coordinates")
g = sns.violinplot(x="Coordinates", y="Unnamed: 0", data=df)
fig = g.get_figure()
fig.savefig("viol.png")
### melt creates long weird df of 2M rows

	



##############################

regScrPrTopicRownames = regScrPrtopic.set_index("Unnamed: 0")
regScrPrTopicRownames.index.names = [None]
regScrPrTopicRownames = regScrPrTopicRownames.drop(["seqnames", "start", "end", "width", "nCounts", "nCells"], axis=1)


regScrPrTopicRownamesFL = regScrPrTopicRownames

for item in regScrList:
	regScrPrTopicRownames = regScrPrTopicRownames.sort_values(by=item, ascending=False)
	regScrPrTopicRownames = regScrPrTopicRownames[[item]] 
	locals()["regScr_"+str(item)] = regScrPrTopicRownames
	regScrPrTopicRownames = regScrPrTopicRownamesFL


#save table as .csv:
regTable.to_csv()

#count neg values:    np.count_nonzero(methTable["DO4341"] < 0)



methTable has: 

                     Unnamed: 0    DO1249    DO1250    DO1253    DO1254    DO1256    DO1257  DO6168    DO6177    DO6195    DO6204           
0          chrY:4868995-4868996  0.040068  0.042150  0.030249  0.045734  0.029226  0.040493  0.060449  0.094056  0.249808  0.097660      
1          chrY:6133739-6133740  0.417496  0.418177  0.399881  0.392938  0.312922  0.431312  0.627849  0.329972  0.432959  0.454886     
2        chrY:22917912-22917913  0.521981  0.494801  0.474009  0.558825  0.520759  0.573331  0.502748  0.591229  0.512314  0.588969     
3          chrY:3446858-3446859  0.874794  0.415923  0.608255  0.435037  0.505973  0.596723  0.922658  0.498599  0.599821  0.421445      
4          chrY:7429348-7429349  0.140750  0.044544  0.068340  0.041774  0.084202  0.032834  0.086345  0.083606  0.090532  0.033507    


Topic_1 = ['DO1398', 'DO1542', 'DO1854', 'DO1879', 'DO1948', 'DO2084', 'DO2096', 'DO2353', 'DO2395', 'DO2860', 'DO3164', 'DO3352', 'DO4149', 'DO44198', 'DO6056']
Topic_2 = ['DO1253', 'DO1262', 'DO1275', 'DO1299', 'DO1342', 'DO1527', 'DO1648', 'DO2180', 'DO2192', 'DO2216', 'DO2234', 'DO2605', 'DO2761', 'DO3025', 'DO3031', 'DO3204', 'DO3626', 'DO4359', 'DO44099', 'DO44133', 'DO44186', 'DO44196', 'DO44200', 'DO4599', 'DO49021', 'DO50021', 'DO5242', 'DO5738', 'DO5759', 'DO6096']



regAssigUnormal has:
                chrY:4868995-4868996  chrY:6133739-6133740  chr22:20342519-20342520  chr22:38598980-38598981  chr22:30112402-30112403
1              	4       		      9                     0                        0                        0
2               0                     8                     0                        0                        0
3               0 	 	              143                   0                        0                        3
4               0                     0                     593                      0                        0


regScrPrtopic has:
                				Scores_Topic3  Scores_Topic4  Scores_Topic5  Scores_Topic6  Scores_Topic7
chrY:4868995-4868996    		0.024327       0.043592       0.006353       0.037375       0.014892	
chrY:6133739-6133740     		0.088911       0.043592       0.006287       0.037375       0.014056
chrY:22917912-22917913     		0.421589       0.043592       0.005740       0.037375       0.019844



topicAssigToPatient has:

									   DO1249    DO1250    DO1253    DO1254    DO1256    DO1257    DO1259
					1  	   			   0.029429  0.008354  0.000003  0.004492  0.032002  0.001979  0.040390
					2  				   0.016377  0.124277  0.018523  0.091671  0.053635  0.062641  0.092560
					3 				   0.185493  0.197464  0.196858  0.191660  0.181314  0.183404  0.178524
					4  				   0.031618  0.000108  0.021195  0.002274  0.028164  0.017091  0.000795
"""

##############################