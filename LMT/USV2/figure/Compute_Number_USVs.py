'''
Created on 12 sept. 2018

@author: Fab
'''


from lmtanalysis.Event import *
from lmtanalysis.Measure import *
import numpy as np; np.random.seed(0)
from tkinter.filedialog import askopenfilename
from lmtanalysis.Util import getMinTMaxTAndFileNameInput, getMinTMaxTInput
import sqlite3
from lmtanalysis.FileUtil import getFilesToProcess, extractPValueFromLMMResult
from lmtanalysis.Animal import AnimalPool
from collections import Counter
from USV.lib.vocUtil import *
import pandas as pd
import seaborn as sns
#import pingouin as pg
from scipy import stats
#from scipy.stats.stats import spearmanr, mannwhitneyu
from statsmodels.stats.anova import AnovaRM
from LMT.USV2.lib.USVUtil import getStrainAgeSexPairGenoPerFile
from LMT.USV2.lib.burster import createBurstFromVoc
from LMT.USV2.figure.figUtil import addJitter, getStarsFromPvalues, strainList,\
    sexList, ageList, acousticVariables
from LMT.USV2.figure.figParameter import colorWT, colorKO, getColorAge,\
    getColorKO, getFigureLabelBurstTraits
'''
from USV.experimentList.experimentList import getExperimentList,\
    getAllExperimentList
'''
from scipy.stats.morestats import wilcoxon
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
#from USV.figure.figUtil import addJitter, getStarsFromPvalues
from USV.figure.figParameter import *
import matplotlib.pyplot as plt
import numpy as np
from USV.burster.burster import *
from scipy.stats import mannwhitneyu, kruskal, ttest_ind
import statsmodels.formula.api as smf
import pandas
'''
from USV.lib.vocUtil import getStrainAgeSexPairGenoPerFile, strainList,\
    sexList, ageList, acousticVariables
'''

measuresList = ['totalNumberUsv', 'totalNumberUsvNight', 'totalNumberUsvDay', 'totalNumberBurst',
                'numberUsvPerBurstList']
burstTrait = [ "durationMs", "nbUSV", "meanDuration", "stdDuration", "meanInterval","stdInterval","meanFrequency" ]
burstTraitList = []
for trait in burstTrait:
    burstTraitList.append('burst'+trait)

measuresList.extend(burstTraitList)


def computeNumberUsv( tmin, tmax, measuresList, strainList, sexList, ageList ):
    
    files = getAllExperimentList()
    #files = getExperimentList(sex='female', age='3mo')
    dataUsv = {}
    for measure in measuresList:
        dataUsv[measure] = {}
        for strain in strainList:
            dataUsv[measure][strain] = {}
            for sex in sexList:
                dataUsv[measure][strain][sex] = {}
                for age in ageList:
                    dataUsv[measure][strain][sex][age] = {}
    
    for file in files:
        expName = file.getFullName()
        print( expName )
        connection = sqlite3.connect( file.file )
        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        #night timeline
        nightTimeLine = EventTimeLine( connection, 'night', idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True )
        print( "loading complete usv timeline dictionary" )
        #usv timeline
        usvTimeLine = EventTimeLine( connection, "Voc", idA=None, minFrame = tmin, maxFrame = tmax, loadEventIndependently=True )
        

        dataUsv['totalNumberUsv'][strainFile][sexFile][ageFile][pairFile] = usvTimeLine.getNumberOfEvent( minFrame = tmin, maxFrame = tmax )

        dataUsv['totalNumberUsvNight'][strainFile][sexFile][ageFile][pairFile] = 0
        for night in nightTimeLine.getEventList():
            nbUsvNight = usvTimeLine.getNumberOfEvent( minFrame=night.startFrame, maxFrame=night.endFrame)
            dataUsv['totalNumberUsvNight'][strainFile][sexFile][ageFile][pairFile] += nbUsvNight

        dataUsv['totalNumberUsvDay'][strainFile][sexFile][ageFile][pairFile] = dataUsv['totalNumberUsv'][strainFile][sexFile][ageFile][pairFile] - dataUsv['totalNumberUsvNight'][strainFile][sexFile][ageFile][pairFile]

        #burst timeline
        burstList = createBurstFromVoc( usvTimeLine )
        dataUsv['totalNumberBurst'][strainFile][sexFile][ageFile][pairFile] = len( burstList )

        for burstTrait in [ "durationMs", "nbUSV", "meanDuration", "stdDuration", "meanInterval","stdInterval","meanFrequency" ]:
            dataUsv['burst'+burstTrait][strainFile][sexFile][ageFile][pairFile] = []

        for burstTrait in ["durationMs", "nbUSV", "meanDuration", "stdDuration", "meanInterval", "stdInterval",
                           "meanFrequency"]:
            for burst in burstList:
                traitForBurst = burst.getValue(burstTrait)
                dataUsv['burst'+burstTrait][strainFile][sexFile][ageFile][pairFile].append( traitForBurst )

        connection.close()

    print(dataUsv.keys())
    with open("dataUsvDescription750.json", 'w') as jFile:
        json.dump(dataUsv, jFile, indent=2)
    print("json file created for USV description")

        
def createDataframeFromJson( jsonFile, strainList, sexList, ageList ):
    #open the json file to work on pre-computed data
    with open( jsonFile ) as json_data:
        dataUsv = json.load(json_data)
    print("json file for acoustic variables re-imported.")
    print(dataUsv.keys())

    constantColumns = [ 'strain', 'sex', 'age', 'pair', 'nbUsv', 'nbUsvDay', 'nbUsvNight', 'nbBurst', 'nbUsvBurst' ]

    dfUsv = pd.DataFrame({}, columns = constantColumns )

    i = 0

    for strain in strainList:
        for sex in sexList:
            for age in ageList:
                for pair in list(dataUsv['totalNumberUsv'][strain][sex][age].keys()):
                    if dataUsv['totalNumberUsv'][strain][sex][age][pair] != None:
                        dfFile = pd.DataFrame(
                            {'strain': strain,
                             'sex': sex,
                             'age': age,
                             'pair': pair,
                             'nbUsv': dataUsv['totalNumberUsv'][strain][sex][age][pair],
                             'nbUsvDay': dataUsv['totalNumberUsvDay'][strain][sex][age][pair],
                             'nbUsvNight': dataUsv['totalNumberUsvNight'][strain][sex][age][pair],
                             'nbBurst': dataUsv['totalNumberBurst'][strain][sex][age][pair],
                             'nbUsvBurst': np.mean(dataUsv['burstnbUSV'][strain][sex][age][pair]),
                            'sdNbUsvBurst': np.std(dataUsv['burstnbUSV'][strain][sex][age][pair])
                            }, index=[i]
                            )
                        dfUsv = pd.concat([dfUsv, dfFile], sort=False)
                        i += 1

    return dfUsv

def createDataframeFromJsonNumberUsvPerBurst( jsonFile, strainList, sexList, ageList ):
    #open the json file to work on pre-computed data
    print( __file__, jsonFile )
    with open( jsonFile ) as json_data:
        print( "caller file: " , __file__ )
        dataUsv = json.load(json_data)
    print("json file for acoustic variables re-imported.")
    print(dataUsv.keys())

    dataDic = { 'strain': [], 'sex': [], 'age': [], 'pair': [], 'nbUsvBurst': [], 'burstmeanDuration': [], 'burstmeanInterval': [] }

    for strain in strainList:
        for sex in sexList:
            for age in ageList:
                for pair in list(dataUsv['totalNumberUsv'][strain][sex][age].keys()):
                    if dataUsv['totalNumberUsv'][strain][sex][age][pair] != None:
                        nbData = len( dataUsv['burstnbUSV'][strain][sex][age][pair] )
                        dataDic['strain'].extend( [strain] * nbData)
                        dataDic['sex'].extend ( [sex] * nbData )
                        dataDic['age'].extend( [age] * nbData )
                        dataDic['pair'].extend( [pair] * nbData )
                        dataDic['nbUsvBurst'].extend( dataUsv['burstnbUSV'][strain][sex][age][pair] )
                        dataDic['burstmeanDuration'].extend(dataUsv['burstmeanDuration'][strain][sex][age][pair])
                        dataDic['burstmeanInterval'].extend(dataUsv['burstmeanInterval'][strain][sex][age][pair])

    df = pd.DataFrame.from_dict(dataDic)
    print('################')
    print(df.head())
    return df

def createDataframeFromJsonAcousticVarAllUsv( jsonFile, acousticVariables, strainList, sexList, ageList ):
    #open the json file to work on pre-computed data
    with open( jsonFile ) as json_data:
        dataUsv = json.load(json_data)
    print("json file for acoustic variables re-imported.")

    colNames = ['strain', 'sex', 'age', 'pair'] + acousticVariables
    print(colNames)
    dataDic = {}
    for col in colNames:
        dataDic[col] = []

    for strain in strainList:
        for sex in sexList:
            for age in ageList:
                for pair in list(dataUsv[strain][sex][age].keys()):
                    for var in acousticVariables:
                        nbData = len( dataUsv[strain][sex][age][pair][var] )
                        dataDic['strain'].extend( [strain] * nbData)
                        dataDic['sex'].extend ( [sex] * nbData )
                        dataDic['age'].extend( [age] * nbData )
                        dataDic['pair'].extend( [pair] * nbData )
                        dataDic[var].extend( dataUsv[strain][sex][age][pair][var] )

    df = pd.DataFrame.from_dict(dataDic)
    print('################')
    print(df.head())
    return df


def getDataFrameWT( dataframe, strain ):
    
    selectedDataframe = dataframe[(dataframe['strain']==strain)]
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print( selectedDataframe.head() )
        
    return selectedDataframe


def getDataFrameKO( dataframe, sex, age ):
    
    selectedDataframe = dataframe[(dataframe['sex']==sex) & (dataframe['age']==age)]
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print( selectedDataframe.head() )
    
    return selectedDataframe


def plotNumberUsvWT(ax, dataframeWT, yMinUsv, yMaxUsv, jitterValue, letter, ageList):
    yLabel = "nb of USVs"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['males', 'females'], rotation=45, FontSize=12, horizontalalignment='right')
    ax.set_ylabel(yLabel, FontSize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv['nbUsv'], yMaxUsv['nbUsv'])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, yMaxUsv['nbUsv'] + 0.06 * (yMaxUsv['nbUsv'] - yMinUsv['nbUsv']), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')

    for age in ageList:
        dataList = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'male')]['nbUsv'].values
        # print(dataList)
        xPositionList = [1] * len(dataList)
        print(dataList)
        ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='v', c=getColorAge(age), alpha=0.5)

        dataList = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'female')]['nbUsv'].values
        # print(dataList)
        xPositionList = [2] * len(dataList)
        print(dataList)
        ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge(age), alpha=0.7)


def plotNumberUsvKO(ax, dataframeKO, yMinUsv, yMaxUsv, jitterValue, letter):
    yLabel = "nb of USVs"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['C57BL/6J', 'Shank3-/-'], rotation=45, fontsize=12, horizontalalignment='right')
    ax.set_ylabel(yLabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in", pad=20)
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv['nbUsv'], yMaxUsv['nbUsv'])
    ax.tick_params(axis='y', labelsize=12)
    #ax.tick_params(axis='x', labelsize=12)
    ax.text(1, -1500, 'F', fontsize=12, horizontalalignment='center')
    ax.text(2, -1500, 'F', fontsize=12, horizontalalignment='center')

    ax.text(-1.2, yMaxUsv['nbUsv'] + 0.06 * (yMaxUsv['nbUsv'] - yMinUsv['nbUsv']), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')

    dataListWT = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'C57BL/6J')]['nbUsv'].values
    # print(dataList)
    xPositionList = [1] * len(dataListWT)
    print(dataListWT)
    ax.scatter(addJitter(xPositionList, jitterValue), dataListWT, marker='o', c=colorWT, alpha=0.9)

    dataListKO = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'Shank3')]['nbUsv'].values
    # print(dataList)
    xPositionList = [2] * len(dataListKO)
    print(dataListKO)
    ax.scatter(addJitter(xPositionList, jitterValue), dataListKO, marker='o', c=colorKO, alpha=0.9)

    U, p = mannwhitneyu(dataListWT, dataListKO, alternative='two-sided')
    print('Nb USVs between genotypes ', ' U = ', U, 'p = ', p)
    ax.text(1.5,
            yMaxUsv['nbUsv'] - 0.06 * (yMaxUsv['nbUsv'] - yMinUsv['nbUsv']), getStarsFromPvalues(p, 1),
            fontsize=20, horizontalalignment='center', color='black')


def plotNumberUsvWTWithAge(ax, dataframeWT, yMinUsv, yMaxUsv, jitterValue, letter, ageList):
    yLabel = "nb of USVs"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1.5, 4.5, 7.5]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['5we', '3mo', '7mo'], rotation=0, fontsize=14, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=14)
    ax.xaxis.set_tick_params(direction="in", pad=20)
    ax.set_xlim(0, 9)
    ax.set_ylim(yMinUsv['nbUsv'], yMaxUsv['nbUsv'])
    ax.tick_params(axis='y', labelsize=12)

    ax.text(-2.5, yMaxUsv['nbUsv'] + 0.06 * (yMaxUsv['nbUsv'] - yMinUsv['nbUsv']), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    xPos = {'5we': {'male':[1], 'female': [2]},
            '3mo': {'male': [4], 'female': [5]},
            '7mo': {'male': [7], 'female': [8]}}
    for age in ageList:
        dataListM = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'male')]['nbUsv'].values
        # print(dataList)
        xPositionList = xPos[age]['male'] * len(dataListM)
        print(dataListM)
        ax.scatter(addJitter(xPositionList, jitterValue), dataListM, marker='v', c=getColorAge(age), alpha=0.9)

        dataListF = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'female')]['nbUsv'].values
        # print(dataList)
        xPositionList = xPos[age]['female'] * len(dataListF)
        print(dataListF)
        ax.scatter(addJitter(xPositionList, jitterValue), dataListF, marker='o', c=getColorAge(age), alpha=0.9)

        ax.text(xPos[age]['male'][0], -1500, 'M', fontsize=12, horizontalalignment='center')
        ax.text(xPos[age]['female'][0], -1500, 'F', fontsize=12, horizontalalignment='center')

        U, p = mannwhitneyu(dataListM, dataListF, alternative='two-sided')
        print('Nb USVs ', age, ' U = ', U, 'p = ', p)
        correction = 3
        if p * correction >= 0.05:
            stars = getStarsFromPvalues(p, 1)
        elif p * correction < 0.05:
            stars = getStarsFromPvalues(p, 1) + '°'
        ax.text((xPos[age]['male'][0]+xPos[age]['female'][0])/2, yMaxUsv['nbUsv'] - 0.06 * (yMaxUsv['nbUsv'] - yMinUsv['nbUsv']), stars,
                    fontsize=20, horizontalalignment='center', color='black', weight='bold')


def plotNumberUsvPerBurstWT (ax, dataframeWT, yMinUsv, yMaxUsv, jitterValue, letter, ageList ):
    yLabel = "nb USV / bursts"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['males', 'females'], rotation=45, FontSize=12, horizontalalignment='right')
    ax.set_ylabel(yLabel, FontSize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv['nbUsvBurst'], yMaxUsv['nbUsvBurst'])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, yMaxUsv['nbUsvBurst'] + 0.06 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')

    for age in ageList:
        dataList = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'male')]['nbUsvBurst'].values
        xPositionList = [1] * len(dataList)
        print(dataList)
        ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='v', c=getColorAge(age),
                   alpha=0.7)

        dataList = dataframeWT[(dataframeWT['age'] == age) & (dataframeWT['sex'] == 'female')]['nbUsvBurst'].values
        xPositionList = [2] * len(dataList)
        print(dataList)
        ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge(age), alpha=0.7)


def plotNumberUsvDayNight( ax, dataframeWT, yMinUsv, yMaxUsv, jitterValue, letter, sexList, ageList ):
    yLabel = "USVs during night (%)"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['C57BL/6J', 'C57BL/6J'], rotation=45, fontsize=12, horizontalalignment='right')
    #ax.set_xticklabels(['M C57BL/6J', 'F C57BL/6J'], rotation=45, FontSize=14, horizontalalignment='right')

    ax.set_ylabel(yLabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in", pad=20)
    ax.set_xlim(0, 2.5)
    #ax.set_ylim(yMinUsv['nbUsv'], yMaxUsv['nbUsv'])
    ax.set_ylim(0, 105)
    ax.tick_params(axis='y', labelsize=12)
    ax.text(1, -8, 'M', fontsize=12, horizontalalignment='center')
    ax.text(2, -8, 'F', fontsize=12, horizontalalignment='center')

    ax.text(-1, 100 + 0.06 * 100, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    xPos = {'male': [1], 'female': [2]}
    marker = {'male': 'v', 'female': 'o'}
    for sex in sexList:
        for age in ageList:
            dataDay = dataframeWT[(dataframeWT['strain'] == 'C57BL/6J') & (dataframeWT['age'] == age) & (dataframeWT['sex'] == sex)]['nbUsvDay'].values
            dataNight = dataframeWT[(dataframeWT['strain'] == 'C57BL/6J') & (dataframeWT['age'] == age) & (dataframeWT['sex'] == sex)]['nbUsvNight'].values
            dataTot = [dataDay[i]+dataNight[i] for i in range(len(dataDay))]
            dataList = [dataNight[i]/dataTot[i]*100 for i in range(len(dataTot))]
            xPositionList = xPos[sex] * len(dataList)
            ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker=marker[sex], c=getColorAge(age), alpha=0.9)

    dataDay = dataframeWT[(dataframeWT['strain'] == 'Shank3') & (dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == 'female')][
        'nbUsvDay'].values
    dataNight = dataframeWT[(dataframeWT['strain'] == 'Shank3') & (dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == 'female')][
        'nbUsvNight'].values
    dataTot = [dataDay[i] + dataNight[i] for i in range(len(dataDay))]
    dataList = [dataNight[i] / dataTot[i] * 100 for i in range(len(dataTot))]
    xPositionList = [3] * len(dataList)
    #ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker=marker['female'], c=getColorKO(), alpha=0.5)


def plotNumberUsvPerBurstKO(ax, dataframeKO, yMinUsv, yMaxUsv, jitterValue, letter ):
    yLabel = "nb USV / bursts"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2]
    ax.set_xticks(xIndex)
    ax.set_xticklabels(['C57BL/6J', 'Shank3'], rotation=45, FontSize=14, horizontalalignment='right')
    ax.set_ylabel(yLabel, FontSize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv['nbUsvBurst'], yMaxUsv['nbUsvBurst'])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, yMaxUsv['nbUsvBurst'] + 0.06 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')


    dataList = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'C57BL/6J')]['nbUsvBurst'].values
    xPositionList = [1] * len(dataList)
    ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge('3mo'), alpha=0.7)

    dataList = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'Shank3')]['nbUsvBurst'].values
    xPositionList = [2] * len(dataList)
    ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorKO(), alpha=0.7)


def plotNumberUsvPerBurstWTBoxplot(ax, dataframe, yMinUsv, yMaxUsv, letter, sexList, ageList):
    yLabel = "nb USV / bursts"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2, 4, 5, 7, 8]
    ax.set_ylabel(yLabel, fontsize=15)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv['nbUsvBurst'], yMaxUsv['nbUsvBurst'])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv['nbUsvBurst'] + 0.06 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='age', y='nbUsvBurst', hue='sex', hue_order=['male', 'female'], data=dataframe, ax=ax, showfliers = False, showmeans=True,
                     notch=True, meanprops=meanprops, width=0.4, dodge=True)
    ax.set_xticklabels(ageList, rotation=0, fontsize=16, horizontalalignment='center')
    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)


    colorList = [getColorAge('5we'), getColorAge('5we'), getColorAge('3mo'), getColorAge('3mo'), getColorAge('7mo'), getColorAge('7mo')]
    edgeList = ['grey', 'black', 'grey', 'black', 'grey', 'black']
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.legend().set_visible(False)

    dataList = {}
    nbBurst = {}
    for sex in sexList:
        dataList[sex] = {}
        nbBurst[sex] = {}
        for age in ageList:
            dataList[sex][age] = dataframe[(dataframe['age'] == age) & (dataframe['sex'] == sex) & (dataframe['strain'] == 'C57BL/6J')]['nbUsvBurst'].values
            nbBurst[sex][age] = len(dataList[sex][age])
            print( sex, age, nbBurst[sex][age], ' bursts')



    xPos = {'5we': 0, '3mo': 1, '7mo': 2}
    for age in ageList:
        '''T, p = ttest_ind(dataList['male'][age], dataList['female'][age])
        print(age, ' male ', nbBurst['male'][age], ' vs female ', nbBurst['female'][age], ': T = ', T, 'p = ', p,
              getStarsFromPvalues(p, 1))'''

        # Mixed model: variable to explain: value; fixed factor = genotype; random effect: group
        dfTestAllAges = dataframe[['age', 'sex', 'pair', 'strain', 'nbUsvBurst']]
        dfTest = dfTestAllAges[dfTestAllAges['age'] == age]
        # create model:
        model = smf.mixedlm("nbUsvBurst ~ sex", dfTest, groups=dfTest['pair'])
        # run model:
        result = model.fit()
        # print summary
        print('Nb USVs per burst', age)
        print(result.summary())
        p, sign = extractPValueFromLMMResult(result=result, keyword='male')
        correction = 3
        if p * correction >= 0.05:
            stars = getStarsFromPvalues(p, 1)
        elif p * correction < 0.05:
            stars = getStarsFromPvalues(p, 1) + '°'

        ax.text(xPos[age],
                yMaxUsv['nbUsvBurst'] - 0.06 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']),
                stars,
                fontsize=16, horizontalalignment='center', color='black', weight='bold')
        ax.text(xPos[age] - 0.1, yMinUsv['nbUsvBurst'] - 0.1 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']), '(n={}) M'.format(nbBurst['male'][age]),
                rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right', color='black')
        ax.text(xPos[age] + 0.2, yMinUsv['nbUsvBurst'] - 0.1 * (yMaxUsv['nbUsvBurst'] - yMinUsv['nbUsvBurst']), '(n={}) F'.format(nbBurst['female'][age]),
                rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right', color='black')


    for sex in sexList:
        T, p = ttest_ind(dataList[sex]['5we'], dataList[sex]['3mo'])
        print(sex, ' 5we ', nbBurst[sex]['5we'], ' vs', sex, ' 3mo ', nbBurst[sex]['3mo'], ': T = ', T, 'p = ', p,
              getStarsFromPvalues(p, 2))
        T, p = ttest_ind(dataList[sex]['3mo'], dataList[sex]['7mo'])
        print(sex, ' 3mo ', nbBurst[sex]['3mo'], ' vs', sex, ' 7mo ', nbBurst[sex]['7mo'], ': T = ', T, 'p = ', p,
              getStarsFromPvalues(p, 2))

        # Mixed model: variable to explain: value; fixed factor = genotype; random effect: group
        dfTestAllAges = dataframe[['age', 'sex', 'pair', 'strain', 'nbUsvBurst']]
        dfTest = dfTestAllAges[dfTestAllAges['sex'] == sex]
        # create model:
        model = smf.mixedlm("nbUsvBurst ~ age", dfTest, groups=dfTest['pair'])
        # run model:
        result = model.fit()
        # print summary
        print('Nb USVs per burst Age effect in ', sex)
        print(result.summary())
        #p, sign = extractPValueFromLMMResult(result=result, keyword='male')






def plotNumberUsvPerBurstKOBoxplot(ax, dataframeKO, burstTrait, yMinUsv, yMaxUsv, letter):
    yLabel = getFigureLabelBurstTraits(burstTrait)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [0, 1]
    ax.set_xticks(xIndex)
    ax.set_ylabel(yLabel, fontsize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 3)
    ax.set_ylim(yMinUsv[burstTrait], yMaxUsv[burstTrait])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, yMaxUsv[burstTrait] + 0.06 * (yMaxUsv[burstTrait] - yMinUsv[burstTrait]), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')
    my_pal = {'C57BL/6J': getColorAge('3mo'), 'Shank3': getColorKO()}
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='strain', y=burstTrait, data=dataframeKO, palette=my_pal, ax=ax, showfliers = False, showmeans=True,
                     notch=True, meanprops=meanprops, width=0.4, dodge=True)
    ax.set_xticklabels(['C57BL/6J', 'Shank3'], rotation=45, fontsize=12, horizontalalignment='right')

    # Add transparency to colors
    for patch in bp.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .7))
    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)

    bp.legend().set_visible(False)

    dataList = {}
    nbBurst = {}

    dataList['B6'] = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'C57BL/6J')][
                burstTrait].values
    dataList['Shank3'] = dataframeKO[(dataframeKO['age'] == '3mo') & (dataframeKO['sex'] == 'female') & (dataframeKO['strain'] == 'Shank3')][
        burstTrait].values
    nbBurst['B6'] = len(dataList['B6'])
    nbBurst['Shank3'] = len(dataList['Shank3'])
    print( 'Shank3: ', nbBurst['Shank3'], 'C57BL/6J: ', nbBurst['B6'])
    print('Shank3: ', np.nanmean(dataList['Shank3']), '+/-', np.nanstd(dataList['Shank3']))
    print('B6: ', np.nanmean(dataList['B6']), '+/-', np.nanstd(dataList['B6']))
    ax.text(0.2, yMinUsv[burstTrait] - 0.1 * (yMaxUsv[burstTrait] - yMinUsv[burstTrait]), '(n={})'.format(nbBurst['B6']),
            rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right', color='black')
    ax.text(1.2, yMinUsv[burstTrait] - 0.1 * (yMaxUsv[burstTrait] - yMinUsv[burstTrait]), '(n={})'.format(nbBurst['Shank3']),
            rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right', color='black')

    '''T, p = ttest_ind(dataList['B6'], dataList['Shank3'])
    print( 'Shank3: ', nbBurst['Shank3'], 'C57BL/6J: ', nbBurst['B6'], ': T = ', T, 'p = ', p, getStarsFromPvalues(p, 1 ))
'''
    dfTest = dataframeKO[['age', 'sex', 'pair', 'strain', burstTrait]]
    dfTest.rename(columns={'age': 'age', 'sex': 'sex', 'pair': 'pair', burstTrait: 'value'}, inplace=True)
    # create model:
    model = smf.mixedlm("value ~ strain", dfTest, groups=dataframeKO['pair'])
    # run model:
    result = model.fit()
    # print summary
    print(result.summary())
    p, sign = extractPValueFromLMMResult(result=result, keyword='Shank3')

    ax.text(0.5,
        yMaxUsv[burstTrait] - 0.06 * (yMaxUsv[burstTrait] - yMinUsv[burstTrait]),
        getStarsFromPvalues(p, 1),
        fontsize=16, horizontalalignment='center', color='black', weight='bold')



if __name__ == '__main__':
    
    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec

    rc('font',**{'family':'serif','serif':['Arial']})
    
    while True:
        
        question = "Do you want to:"
        question +="\n\t [c]ompute the total number of USVs, the number of USV bursts and the number of USV per burst over the three days?"
        question +="\n\t [p]lot figure 2 with data for C57BL/6J mice and for Shank3 mice?"
        question +="\n\t [e]xtract the mean and sd for different categories?"
        question +="\n"
        answer = input("Action:" )
        
        if answer=="c":
                        
            computeNumberUsv( tmin=0, tmax=7776000, measureList=measuresList, strainList=strainList, sexList=sexList, ageList=ageList )
            print("ok")
            
            break

        
        if answer=="p":
            
            df = createDataframeFromJson( jsonFile = 'dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList )
            #print(df)
            ###################################################
            #data wild-type
            dataframeWT = getDataFrameWT( df, strain='C57BL/6J' )

            #################################################
            for var in ['nbUsv', 'nbUsvDay', 'nbUsvNight', 'nbBurst', 'nbUsvBurst']:
                print('tests for age effect: ')
                for sexClass in sexList:
                    K, p = stats.kruskal(dataframeWT[var][(dataframeWT['age']=='5we') & (dataframeWT['sex']==sexClass)],
                                  dataframeWT[var][(dataframeWT['age']=='3mo') & (dataframeWT['sex']==sexClass)],
                                  dataframeWT[var][(dataframeWT['age']=='7mo') & (dataframeWT['sex']==sexClass)])
                    print(var, ' ', sexClass, ' K = ', K, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age']=='5we') & (dataframeWT['sex']==sexClass)], dataframeWT[var][(dataframeWT['age']=='3mo') & (dataframeWT['sex']==sexClass)], alternative='two-sided')
                    print(var, ' ', sexClass, ' 5we vs 3mo', ' U = ', U, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age']=='3mo') & (dataframeWT['sex']==sexClass)], dataframeWT[var][(dataframeWT['age']=='7mo') & (dataframeWT['sex']==sexClass)], alternative='two-sided')
                    print(var, ' ', sexClass, ' 3mo vs 7mo', ' U = ', U, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age']=='5we') & (dataframeWT['sex']==sexClass)], dataframeWT[var][(dataframeWT['age']=='7mo') & (dataframeWT['sex']==sexClass)], alternative='two-sided')
                    print(var, ' ', sexClass, ' 5we vs 7mo', ' U = ', U, 'p = ', p)
                
                print('Tests between sexes:')
                for ageClass in ageList:
                    U, p = mannwhitneyu(dataframeWT[var][(dataframeWT['age']==ageClass) & (dataframeWT['sex']=='male')], dataframeWT[var][(dataframeWT['age']==ageClass) & (dataframeWT['sex']=='female')], alternative='two-sided')
                    print(var, ' ', ageClass, ' U = ', U, 'p = ', p)
                
            #########################################################
            #plots
            yMinUsv = {'nbUsv': 0, 'nbBurst': 0, 'nbUsvBurst': 0}
            yMaxUsv = {'nbUsv': 24000, 'nbBurst': 5000, 'nbUsvBurst': 50}
            jitterValue = 0.2
            x1 = {'nbUsvDay': {'5we': 1, '3mo': 4, '7mo': 7}, 'nbUsvNight': {'5we': 2, '3mo': 5, '7mo': 8},
                  'nbBurst': {'5we': 1.5, '3mo': 4.5, '7mo': 7.5}, 'nbUsvBurst':{'5we': 1.5, '3mo': 4.5, '7mo': 7.5} }

            # Create general figure
            gs = gridspec.GridSpec(1, 7)
            # fig = plt.subplots(figsize=(24, 15), sharex=False, sharey=False)
            fig = plt.figure(figsize=(14, 3))

            ################################################
            # plot the number of usv in males versus females with all ages
            #plotNumberUsvWT( ax = fig.add_subplot(gs[0, 0:1]), dataframeWT = dataframeWT, yMinUsv = yMinUsv, yMaxUsv = yMaxUsv, jitterValue = jitterValue, letter = 'A' )
            plotNumberUsvWTWithAge(ax=fig.add_subplot(gs[0, 0:2]), dataframeWT=dataframeWT, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv,
                            jitterValue=jitterValue, letter='A', ageList=ageList)

            ############################################################
            # plot the number of usv within bursts for males and females
            #plotNumberUsvPerBurstWT( ax = fig.add_subplot(gs[0, 1:2]), dataframeWT = dataframeWT, yMinUsv = yMinUsv, yMaxUsv = yMaxUsv, jitterValue=jitterValue, letter = 'B )
            # plot the number of USV within bursts as boxplots for males and females across age classes
            df = createDataframeFromJsonNumberUsvPerBurst(jsonFile='dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList)
            dataframeWT = getDataFrameWT(df, strain = 'C57BL/6J')

            plotNumberUsvPerBurstWTBoxplot(ax=fig.add_subplot(gs[0, 2:4]), dataframe=dataframeWT, yMinUsv=yMinUsv,
                                           yMaxUsv=yMaxUsv, letter = 'B', sexList=sexList, ageList=ageList)
            print('boxplot1 done.')

            ################################################
            # plot the number of usv in days and nights
            df = createDataframeFromJson( jsonFile = 'dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList )
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            plotNumberUsvDayNight( ax = fig.add_subplot(gs[0, 4:5]), dataframeWT = dataframeWT, yMinUsv = yMinUsv, yMaxUsv = yMaxUsv, jitterValue=jitterValue, letter = 'D', sexList=sexList, ageList=ageList )
            print('day night done')
            ###################################################
            # data KO with WT controls
            dataframeKO = getDataFrameKO(df, sex = 'female', age = '3mo')
            # run non parametric test:
            dfKO = dataframeKO[dataframeKO['strain'] == 'Shank3']
            dfWT = dataframeKO[dataframeKO['strain'] == 'C57BL/6J']
            for var in ['nbUsv', 'nbBurst', 'nbUsvBurst']:
                U, p = mannwhitneyu(dfKO[var], dfWT[var], alternative='two-sided')
                print(var, ' U = ', U, 'p = ', p)
            '''
            # plot the number of usv within bursts for WT and KO females
            plotNumberUsvPerBurstKO(ax=fig.add_subplot(gs[0, 5:6]), dataframeKO=dataframeKO, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, jitterValue=jitterValue, letter = 'E' )
            '''
            #plot the number of USV within bursts as boxplots
            df = createDataframeFromJsonNumberUsvPerBurst( jsonFile = 'dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList )
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')

            plotNumberUsvPerBurstKOBoxplot(ax=fig.add_subplot(gs[0, 5:6]), dataframeKO=dataframeKO, burstTrait='nbUsvBurst', yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, letter = 'E')



            plt.tight_layout()
            plt.show()
            print("Job done.")

            break

        if answer == "e":

            df = createDataframeFromJson(jsonFile='dataUsvDescription750.json')
            print(df)

            for var in ['nbUsv', 'nbUsvDay', 'nbUsvNight', 'nbBurst', 'nbUsvBurst']:
                dataWt = df[var][(df['strain'] == 'C57BL/6J') & (df['age'] == '3mo') & (df['sex'] == 'female')]
                dataKo = df[var][(df['strain'] == 'Shank3') & (df['age'] == '3mo') & (df['sex'] == 'female')]

                meanWt = np.mean(dataWt)
                meanKo = np.mean(dataKo)
                sdWt = np.std(dataWt)
                sdKo = np.std(dataKo)
                print('{} WT: {} +/- {} KO: : {} +/- {}'.format( var, meanWt, sdWt, meanKo, sdKo ) )

            quit()
            ###################################################
            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            #################################################
            for var in ['nbUsv', 'nbUsvDay', 'nbUsvNight', 'nbBurst', 'nbUsvBurst']:
                print('tests for age effect: ')
                for sexClass in sexList:
                    K, p = stats.kruskal(
                        dataframeWT[var][(dataframeWT['age'] == '5we') & (dataframeWT['sex'] == sexClass)],
                        dataframeWT[var][(dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == sexClass)],
                        dataframeWT[var][(dataframeWT['age'] == '7mo') & (dataframeWT['sex'] == sexClass)])

            break

        if answer == "pv":
            df = createDataframeFromJsonAcousticVarAllUsv(jsonFile='dataAcousticAnalysisAllUsvs_all_pairs.json', acousticVariables=acousticVariables, strainList=strainList, sexList=sexList, ageList=ageList)
            print(df.head())
            break