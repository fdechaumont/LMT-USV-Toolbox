'''
Created on 09 Dec. 2020

@author: Elodie
'''

from lmtanalysis.Event import *
from lmtanalysis.Measure import *
import numpy as np; np.random.seed(0)
from tkinter.filedialog import askopenfilename
from lmtanalysis.Util import getMinTMaxTAndFileNameInput, getMinTMaxTInput
import sqlite3
from lmtanalysis.FileUtil import getFilesToProcess, getFigureBehaviouralEventsLabelsFrench
from lmtanalysis.Animal import AnimalPool
from collections import Counter
from USV.lib.vocUtil import *
import pandas as pd
import seaborn as sns
#import pingouin as pg
from scipy import stats
#from scipy.stats.stats import spearmanr, mannwhitneyu
from statsmodels.stats.anova import AnovaRM

from scipy.stats.morestats import wilcoxon
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

from USV.figure.figParameter import *
import matplotlib.pyplot as plt
import numpy as np
from USV.burster.burster import *
#from USV.usvDescription.Compute_Number_USVs import getDataFrameWT
from USV.figure.figUtil import *
from LMT.USV2.lib.USVUtil import getStrainAgeSexPairGenoPerFile
from LMT.USV2.lib.burster import createBurstFromVoc
from LMT.USV2.figure.figParameter import getFigureBehaviouralEvents,\
    getFigureBehaviouralEventsLabels, getColorAge, getColorKO, getColorGeno
from LMT.USV2.figure.figUtil import strainList, sexList, ageList, addJitter,\
    getStarsFromPvalues
from scipy.stats.mstats_basic import mannwhitneyu
from LMT.USV2.experimentList.experimentList import getAllExperimentList
from LMT.USV2.figure.Compute_Number_USVs import getDataFrameWT


def eventMatchWithTimeLine( eventATimeLine, eventBTimeLine):
    # This function checks how many eventB overlap with eventA
    ratio = eventMatchWithDic( eventATimeLine.getEventList(), eventBTimeLine.getEventList() )
    return ratio


def eventMatchWithDic( eventListA, eventListB):
    #This function checks how many eventB overlap with eventA
    #create dictionary from the list of event A (the event to be checked)
    dicA = {}
    for event in eventListA:
        try:
            for t in range( event.startFrame, event.endFrame+1):
                dicA[t] = True
        except:
            for t in range( event.getStartFrame(), event.getEndFrame()+1):
                dicA[t] = True
    #check for correlation with the list of event B
    nbMatch = 0

    for event in eventListB:
        try:
            for t in range( event.startFrame, event.endFrame+1):
                if t in dicA.keys():
                    nbMatch += 1
                    break
        except:
            for t in range( event.getStartFrame(), event.getEndFrame()+1):
                if t in dicA.keys():
                    nbMatch += 1
                    break

    ratio = nbMatch / len( eventListB)
    return ratio


def computeCorrelationBurstWithEvent( files, tmin, tmax, eventListToTest, strainList, sexList, ageList ):
    '''
    This code focuses on the USV / bursts and investigates with which events they are overlapping.
    This analysis used the USV / bursts as the reference events.
    It builds a json file to store data.
    '''
    data = {}
    for cat in ['usv', 'burst']:
        data[cat] = {}
        for event in eventListToTest:
            data[cat][event] = {}
            for strain in strainList:
                data[cat][event][strain] = {}
                for sex in sexList:
                    data[cat][event][strain][sex] = {}
                    for age in ageList:
                        data[cat][event][strain][sex][age] = {}

    for file in files:
        expName = file.getFullName()
        print( expName )
        connection = sqlite3.connect( file.file )
        pool = AnimalPool()
        pool.loadAnimals(connection)

        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        #generate usv timeline
        usvTimeLine = EventTimeLine(connection, "Voc", idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        
        #create burst list based on USV
        burstList = createBurstFromVoc(usvTimeLine)

        for eventToTest in eventListToTest:
            eventToTestTimeLine = {}
            data['burst'][eventToTest][strainFile][sexFile][ageFile][pairFile] = {}
            data['usv'][eventToTest][strainFile][sexFile][ageFile][pairFile] = {}

            for animal in pool.getAnimalList():
                # Load the timeline for the event to test for each animal
                print("loading  for animal: ", animal.RFID)
                eventToTestTimeLine[animal] = EventTimeLine(connection, eventToTest, idA=animal.baseId, minFrame=tmin, maxFrame=tmax)
                ratioBurst = eventMatchWithDic( eventToTestTimeLine[animal].eventList, burstList )
                ratioUsv = eventMatchWithTimeLine( eventToTestTimeLine[animal], usvTimeLine )

                data['burst'][eventToTest][strainFile][sexFile][ageFile][pairFile][animal.RFID] = ratioBurst
                data['usv'][eventToTest][strainFile][sexFile][ageFile][pairFile][animal.RFID] = ratioUsv
                print(eventToTest, strainFile, sexFile, ageFile, pairFile, animal.RFID, ratioBurst, ratioUsv)
            print('#####################################')
        connection.close()

    return data

def computeCorrelationBurstWithEventDiffGeno( files, tmin, tmax, eventListToTest, strainList, sexList, genoList ):
    '''
    This code focuses on the USV / bursts and investigates with which events they are overlapping.
    This analysis used the USV / bursts as the reference events.
    It builds a json file to store data.
    '''
    data = {}
    for cat in ['usv', 'burst']:
        data[cat] = {}
        for event in eventListToTest:
            data[cat][event] = {}
            for strain in strainList:
                data[cat][event][strain] = {}
                for sex in sexList:
                    data[cat][event][strain][sex] = {}
                    for geno in genoList:
                        data[cat][event][strain][sex][geno] = {}


    for file in files:

        connection = sqlite3.connect( file )
        pool = AnimalPool()
        pool.loadAnimals(connection)

        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        #generate usv timeline
        usvTimeLine = EventTimeLine(connection, "Voc", idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        #cleanVoc(usvTimeLine)
        #create burst list based on USV
        burstList = createBurstFromVoc(usvTimeLine)

        for eventToTest in eventListToTest:
            eventToTestTimeLine = {}
            data['burst'][eventToTest][strainFile][sexFile][genoFile][pairFile] = {}
            data['usv'][eventToTest][strainFile][sexFile][genoFile][pairFile] = {}

            for animal in pool.getAnimalList():
                # Load the timeline for the event to test for each animal
                print("loading  for animal: ", animal.RFID)
                eventToTestTimeLine[animal] = EventTimeLine(connection, eventToTest, idA=animal.baseId, minFrame=tmin, maxFrame=tmax)
                ratioBurst = eventMatchWithDic( eventToTestTimeLine[animal].eventList, burstList )
                ratioUsv = eventMatchWithTimeLine( eventToTestTimeLine[animal], usvTimeLine )

                data['burst'][eventToTest][strainFile][sexFile][genoFile][pairFile][animal.RFID] = ratioBurst
                data['usv'][eventToTest][strainFile][sexFile][genoFile][pairFile][animal.RFID] = ratioUsv
                print(eventToTest, strainFile, sexFile, genoFile, pairFile, animal.RFID, ratioBurst, ratioUsv)
            print('#####################################')

        connection.close()
    return data


def computeCorrelationEventWithBurst( files, tmin, tmax, eventListToTest, strainList, sexList, ageList ):
    '''
    This code focuses on the events and investigates with which USVs / bursts they are overlapping.
    This analysis used the events as the reference events.
    It builds a json file to store data.
    '''
    data = {}
    for cat in ['usv', 'burst']:
        data[cat] = {}
        for event in eventListToTest:
            data[cat][event] = {}
            for strain in strainList:
                data[cat][event][strain] = {}
                for sex in sexList:
                    data[cat][event][strain][sex] = {}
                    for age in ageList:
                        data[cat][event][strain][sex][age] = {}

    for file in files:
        expName = file.getFullName()
        print( expName )
        connection = sqlite3.connect( file.file )
        pool = AnimalPool()
        pool.loadAnimals(connection)

        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        #generate usv timeline
        usvTimeLine = EventTimeLine(connection, "Voc", idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        
        #create burst list based on USV
        burstList = createBurstFromVoc(usvTimeLine)

        for eventToTest in eventListToTest:
            eventToTestTimeLine = {}
            data['burst'][eventToTest][strainFile][sexFile][ageFile][pairFile] = {}
            data['usv'][eventToTest][strainFile][sexFile][ageFile][pairFile] = {}

            for animal in pool.getAnimalList():
                # Load the timeline for the event to test for each animal
                print("loading  for animal: ", animal.RFID)
                eventToTestTimeLine[animal] = EventTimeLine(connection, eventToTest, idA=animal.baseId, minFrame=tmin, maxFrame=tmax)
                ratioBurst = eventMatchWithDic( burstList, eventToTestTimeLine[animal].eventList )
                ratioUsv = eventMatchWithTimeLine( usvTimeLine, eventToTestTimeLine[animal] )

                data['burst'][eventToTest][strainFile][sexFile][ageFile][pairFile][animal.RFID] = ratioBurst
                data['usv'][eventToTest][strainFile][sexFile][ageFile][pairFile][animal.RFID] = ratioUsv
                print(eventToTest, strainFile, sexFile, ageFile, pairFile, animal.RFID, ratioBurst, ratioUsv)
            print('#####################################')

        connection.close()
    return data

def computeCorrelationEventWithBurstDiffGeno( files, tmin, tmax, eventListToTest, strainList, sexList, genoList ):
    '''
    This code focuses on the events and investigates with which USVs / bursts they are overlapping.
    This analysis used the events as the reference events.
    It builds a json file to store data.
    '''
    data = {}
    for cat in ['usv', 'burst']:
        data[cat] = {}
        for event in eventListToTest:
            data[cat][event] = {}
            for strain in strainList:
                data[cat][event][strain] = {}
                for sex in sexList:
                    data[cat][event][strain][sex] = {}
                    for geno in genoList:
                        data[cat][event][strain][sex][geno] = {}

    for file in files:
        connection = sqlite3.connect( file )
        pool = AnimalPool()
        pool.loadAnimals(connection)

        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        #generate usv timeline
        usvTimeLine = EventTimeLine(connection, "Voc", idA=None, minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        #cleanVoc(usvTimeLine)
        #create burst list based on USV
        burstList = createBurstFromVoc(usvTimeLine)

        for eventToTest in eventListToTest:
            eventToTestTimeLine = {}
            data['burst'][eventToTest][strainFile][sexFile][genoFile][pairFile] = {}
            data['usv'][eventToTest][strainFile][sexFile][genoFile][pairFile] = {}

            for animal in pool.getAnimalList():
                # Load the timeline for the event to test for each animal
                print("loading  for animal: ", animal.RFID)
                eventToTestTimeLine[animal] = EventTimeLine(connection, eventToTest, idA=animal.baseId, minFrame=tmin, maxFrame=tmax)
                ratioBurst = eventMatchWithDic( burstList, eventToTestTimeLine[animal].eventList )
                ratioUsv = eventMatchWithTimeLine( usvTimeLine, eventToTestTimeLine[animal] )

                data['burst'][eventToTest][strainFile][sexFile][genoFile][pairFile][animal.RFID] = ratioBurst
                data['usv'][eventToTest][strainFile][sexFile][genoFile][pairFile][animal.RFID] = ratioUsv
                print(eventToTest, strainFile, sexFile, genoFile, pairFile, animal.RFID, ratioBurst, ratioUsv)
            print('#####################################')

        connection.close()
    return data


def saveResultsAsJson( jsonFile, resultFile):
    with open(jsonFile, 'w') as jFile:
        json.dump(resultFile, jFile, indent=4)
    print("json file created")


def generateDataframeCorrelationFromDic( data ):
    #generate data frame from dictionary:
    constantColumns = [ 'strain', 'sex', 'age', 'pair', 'event', 'usvCorr', 'burstCorr' ]
    eventListToTest = getFigureBehaviouralEvents(longList=True)
    df = pd.DataFrame({}, columns = constantColumns )

    i = 0
    for event in eventListToTest:
        for strain in strainList:
            for sex in sexList:
                for age in ageList:
                    for pair in list(data['burst'][event][strain][sex][age].keys()):
                        for animal in list(data['burst'][event][strain][sex][age][pair].keys()):
                            #if (data['burst'][event][strain][sex][age][pair][animal] != None) & (data['usv'][event][strain][sex][age][pair][animal] != None):
                            dfFile = pd.DataFrame(
                                {'strain': strain,
                                 'sex': sex,
                                 'age': age,
                                 'pair': pair,
                                 'event': event,
                                 'usvCorr': data['usv'][event][strain][sex][age][pair][animal],
                                 'burstCorr': data['burst'][event][strain][sex][age][pair][animal]
                                }, index=[i]
                                )
                            df = pd.concat([df, dfFile], sort=False)
                            i += 1

    print(df.head())
    return df

def generateDataframeCorrelationFromDicDiffGeno( data, eventListToTest, strainList, sexList, genoList ):
    #generate data frame from dictionary:
    constantColumns = [ 'strain', 'sex', 'geno', 'pair', 'event', 'usvCorr', 'burstCorr' ]
    df = pd.DataFrame({}, columns = constantColumns )

    i = 0
    for event in eventListToTest:
        for strain in strainList:
            for sex in sexList:
                for geno in genoList:
                    for pair in list(data['burst'][event][strain][sex][geno].keys()):
                        for animal in list(data['burst'][event][strain][sex][geno][pair].keys()):
                            #if (data['burst'][event][strain][sex][age][pair][animal] != None) & (data['usv'][event][strain][sex][age][pair][animal] != None):
                            dfFile = pd.DataFrame(
                                {'strain': strain,
                                 'sex': sex,
                                 'geno': geno,
                                 'pair': pair,
                                 'event': event,
                                 'usvCorr': data['usv'][event][strain][sex][geno][pair][animal],
                                 'burstCorr': data['burst'][event][strain][sex][geno][pair][animal]
                                }, index=[i]
                                )
                            df = pd.concat([df, dfFile], sort=False)
                            i += 1

    print(df.head())
    return df

def plotCorrelationUsvEvents(ax, selectedEventList, dataframe, cat, yMin, yMax, ageList, jitterValue, letter ):
    yLabel = '% USVs with events'
    eventLabelList = []
    for event in selectedEventList:
        eventLabelList.append( getFigureBehaviouralEventsLabels(event))

    xPos = [2, 7, 12]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = xPos
    ax.set_xticks(xIndex)
    ax.set_xticklabels(eventLabelList, rotation=45, fontsize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 14)
    ax.set_ylim(yMin, yMax)
    ax.tick_params(axis='y', labelsize=12)

    ax.text(-2.5, yMax + 0.06 * (yMax - yMin), letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    jitterValue = 0.3
    n = 0
    for event in selectedEventList:
        for age in ageList:
            dataList = dataframe[(dataframe['strain'] == 'C57BL/6J') & (dataframe['event'] == event) & (dataframe['age'] == age) & (dataframe['sex'] == 'male')][cat].values
            dataList = [dataList[i]*100 for i in range(len(dataList))]
            x = xPos[n]-1
            xPositionList = [x] * len(dataList)
            print(dataList)
            ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='v', c=getColorAge(age), alpha=0.7)

            dataList = dataframe[(dataframe['strain'] == 'C57BL/6J') & (dataframe['event'] == event) & (dataframe['age'] == age) & (dataframe['sex'] == 'female')][cat].values
            dataList = [dataList[i]*100 for i in range(len(dataList))]
            x = xPos[n]
            xPositionList = [x] * len(dataList)
            print(dataList)
            ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge(age), alpha=0.7)

        dataList = dataframe[(dataframe['strain'] == 'Shank3') & (dataframe['event'] == event) & (dataframe['age'] == '3mo') & (dataframe['sex'] == 'female')][cat].values
        dataList = [dataList[i] * 100 for i in range(len(dataList))]
        x = xPos[n]+1
        xPositionList = [x] * len(dataList)
        print(dataList)
        ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorKO(), alpha=0.7)

        n += 1


def plotCorrelationUsvEventsFemalesOnly(ax, selectedEventList, dataframe, cat, yMin, yMax, ageList, jitterValue, letter ):
    yLabel = "proportion {}".format(cat)
    eventLabelList = []
    for event in selectedEventList:
        eventLabelList.append(getFigureBehaviouralEventsLabels(event))
    xPos = [1, 3, 5, 7]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = xPos
    ax.set_xticks(xIndex)
    ax.set_xticklabels(eventLabelList, rotation=45, FontSize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, FontSize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 8)
    ax.set_ylim(yMin, yMax)
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2, yMax + 0.06 * (yMax - yMin), letter, FontSize=20, horizontalalignment='center', color='black', weight='bold')

    n = 0
    for event in selectedEventList:
        for age in ageList:
            dataList = dataframe[(dataframe['event'] == event) & (dataframe['age'] == age) & (dataframe['sex'] == 'female')][cat].values
            x = xPos[n]
            xPositionList = [x] * len(dataList)
            print(dataList)
            ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge(age), alpha=0.7)
        n += 1

def plotCorrelationUsvEventsFemalesOnlyDiffGeno(ax, selectedEventList, dataframe, cat, yMin, yMax, genoList, jitterValue, letter, ylabel ):
    yLabel = ylabel
    eventLabelList = []
    for event in selectedEventList:
        eventLabelList.append(getFigureBehaviouralEventsLabelsFrench(event))
    xPos = [1, 3, 5, 7, 9, 11, 13, 15, 17]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = xPos
    ax.set_xticks(xIndex)
    ax.set_xticklabels(eventLabelList, rotation=45, fontsize=12, horizontalalignment='right')
    ax.set_ylabel(yLabel, fontsize=15)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 18)
    ax.set_ylim(yMin, yMax)
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-2.5, yMax + 0.06 * (yMax - yMin), letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')

    d = {genoList[0]: -0.3, genoList[1]: 0.3}
    n = 0
    for event in selectedEventList:
        dataDic = {}
        for geno in genoList:
            dataDic[geno] = dataframe[(dataframe['event'] == event) & (dataframe['geno'] == geno) & (dataframe['sex'] == 'female')][cat].values
            x = xPos[n] + d[geno]
            xPositionList = [x] * len(dataDic[geno])
            print(dataDic[geno])
            ax.scatter(addJitter(xPositionList, jitterValue), dataDic[geno], marker='o', c=getColorGeno(geno), alpha=0.7)

        ax.errorbar(xPos[n] + d[genoList[0]], np.mean(dataDic[genoList[0]]), yerr=np.std(dataDic[genoList[0]]), marker='o', c='black')
        ax.errorbar(xPos[n] + d[genoList[1]], np.mean(dataDic[genoList[1]]), yerr=np.std(dataDic[genoList[1]]), marker='o', c='black')
        W, p = mannwhitneyu( dataDic[genoList[0]], dataDic[genoList[1]])
        print( 'Mann-Whitney U for correlation USV vs {}, {} n={}, {} n={}: W={}, p={}'.format( event, genoList[0], len(dataDic[genoList[0]]), genoList[1], len(dataDic[genoList[1]]), W, p ))
        ax.text(xPos[n], 0.7, getStarsFromPvalues(pvalue=p, numberOfTests=len(selectedEventList)), fontsize=12)

        n += 1



def plotCorrelationUsvEventsWithUSVWithKo(ax, selectedEventList, dataframe, cat, yMin, yMax, ageList, jitterValue, letter ):
    yLabel = '% events with USVs'
    eventLabelList = []
    for event in selectedEventList:
        eventLabelList.append(getFigureBehaviouralEventsLabels(event))

    xPos = [2, 7, 12]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = xPos
    ax.set_xticks(xIndex)
    #ax.set_xticklabels(selectedEventList, rotation=45, FontSize=12, horizontalalignment='center')
    ax.set_xticklabels(eventLabelList, rotation=45, fontsize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 14)
    ax.set_ylim(yMin, yMax)
    ax.tick_params(axis='y', labelsize=12)

    ax.text(-3, yMax + 0.06 * (yMax - yMin), letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')

    n = 0
    for event in selectedEventList:
        for age in ageList:
            dataList = dataframe[(dataframe['event'] == event) & (dataframe['age'] == age) & (dataframe['sex'] == 'female') & (dataframe['strain'] == 'C57BL/6J') ][cat].values
            dataList = [dataList[i]*100 for i in range(len(dataList))]
            x = xPos[n]-1
            xPositionList = [x] * len(dataList)
            print(dataList)
            ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorAge(age), alpha=0.7)
            if age == '3mo':
                dataList = dataframe[(dataframe['event'] == event) & (dataframe['age'] == age) & (dataframe['sex'] == 'female') & (dataframe['strain'] == 'Shank3') ][cat].values
                dataList = [dataList[i] * 100 for i in range(len(dataList))]
                x = xPos[n] + 1
                xPositionList = [x] * len(dataList)
                print(dataList)
                ax.scatter(addJitter(xPositionList, jitterValue), dataList, marker='o', c=getColorKO(), alpha=0.7)
        n += 1


if __name__ == '__main__':
    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})

    eventListToTest = getFigureBehaviouralEvents(longList=True)

    while True:

        question = "Do you want to:"
        question += "\n\t [cv]ompute the correlations between USV/bursts and events?"
        question += "\n\t [ce]ompute the correlations between events and USV/bursts?"
        question += "\n\t [p]lot data?"
        question += "\n"
        answer = input(question)

        if answer == "cv":
            jsonFile = "dataCorrelationRefUsvBurst.json"
            files = getAllExperimentList()
            tmin = 0
            tmax = 7776000
            data = computeCorrelationBurstWithEvent( files = files, tmin = tmin, tmax = tmax, eventListToTest = eventListToTest, strainList = strainList, sexList = sexList, ageList = ageList )

            saveResultsAsJson( jsonFile, data )

            break

        if answer == "ce":
            jsonFile = "dataCorrelationRefEvent.json"
            files = getAllExperimentList()
            tmin = 0
            tmax = 7776000
            data = computeCorrelationEventWithBurst( files = files, tmin = tmin, tmax = tmax, eventListToTest = eventListToTest, strainList = strainList, sexList = sexList, ageList = ageList )

            saveResultsAsJson( jsonFile, data )

            break

        if answer == "p":
            selectedEventList = ['Stop isolated', 'Contact', 'Train2', 'longChase']
            # Create general figure
            gs = gridspec.GridSpec(1, 8)
            jitterValue = 0.7
            fig = plt.figure(figsize=(16, 3))

            #plot correlation of USV / bursts with selected behavioural events
            # open the json file to work on pre-computed data
            jsonFile = "dataCorrelationRefUsvBurst_2000.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file for acoustic variables re-imported.")
            #generate data frame from dictionary:
            df = generateDataframeCorrelationFromDic( data )

            df = getDataFrameWT(df, 'C57BL/6J')

            plotCorrelationUsvEvents(ax = fig.add_subplot(gs[0, 0:2]), selectedEventList = selectedEventList, dataframe = df, cat = 'usvCorr', yMin=0, yMax=1, ageList=ageList, jitterValue = jitterValue,
                                     letter = 'A')

            plotCorrelationUsvEvents(ax=fig.add_subplot(gs[0, 2:4]),
                                     selectedEventList=selectedEventList,
                                     dataframe=df, cat='burstCorr', yMin=0, yMax=1, ageList=ageList,
                                     jitterValue=jitterValue,
                                     letter='B')

            # plot correlation of selected behavioural events with USV / bursts
            selectedEventList = ['Oral-genital Contact', 'FollowZone Isolated', 'Train2', 'longChase']
            # open the json file to work on pre-computed data
            jsonFile = "dataCorrelationRefEvent_2000.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file re-imported.")
            # generate data frame from dictionary:
            df = generateDataframeCorrelationFromDic(data)

            df = getDataFrameWT(df, 'C57BL/6J')

            plotCorrelationUsvEventsFemalesOnly(ax=fig.add_subplot(gs[0, 4:6]), selectedEventList=selectedEventList, dataframe=df,
                                     cat='usvCorr', yMin=0, yMax=1, ageList=ageList, jitterValue=jitterValue,
                                     letter='C')

            plotCorrelationUsvEventsFemalesOnly(ax=fig.add_subplot(gs[0, 6:8]),
                                     selectedEventList=selectedEventList,
                                     dataframe=df, cat='burstCorr', yMin=0, yMax=1, ageList=ageList,
                                     jitterValue=jitterValue,
                                     letter='D')
            plt.tight_layout()
            plt.show()
            print("Job done.")


            break


