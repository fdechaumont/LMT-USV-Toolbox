'''
Created on 20 janv. 2020

@author: Elodie Ey
'''




from tkinter.filedialog import askopenfilename
import sqlite3
import os
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

#from USV.figure.featureHeatMap import featureHeatMapPValLegend, featureHeatMapEffectSizeLegend
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Event import Event, EventTimeLine
from lmtanalysis.Measure import oneDay, oneHour, oneMinute
from lmtanalysis.Util import getMinTMaxTAndFileNameInput, convert_to_d_h_m_s,\
    getMinTMaxTInput
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
import pandas
import scipy.stats as stats
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu


import numpy as np

from scipy import signal
from scipy.io import wavfile
import os
import wave
import pylab
from lmtanalysis.FileUtil import getFilesToProcess, extractPValueFromLMMResult
import pandas as pd
from pandas.core.frame import DataFrame
from scipy import stats

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import json
from LMT.USV2.figure.figUtil import strainList, sexList, yMaxVar, yMinVar,\
    getStarsFromPvalues
from LMT.USV2.figure.Compute_Speed_Duration_Events_With_USV import ageList
from LMT.USV2.lib.USVUtil import getStrainAgeSexPairGenoPerFile
from LMT.USV2.figure.figParameter import getFigureLabelTrait,\
    getFigureBehaviouralEventsLabels, getColorAge, getColorWT, getColorKO,\
    getColorEvent, getFigureBehaviouralEvents, getFigureVocTraits
from LMT.USV2.figure.Compute_Number_USVs import getDataFrameWT, getDataFrameKO

'''
from LMT.USV.figure.vocTraitUsagePerEventContext import plotVocTraitUsagePerEventContextPerSetShort2
from LMT.USV.lib.vocUtil import *
from LMT.USV.experimentList.experimentList import getAllExperimentList, getExperimentList
from LMT.USV.figure.figParameter import getFigureBehaviouralEvents, colorWT,\
    colorKO, getPaperColor, getFigureVocTraits, getFigureLabelTrait,\
    getFigureBehaviouralEventsLabels, getColorAge, getColorWT, getColorKO,\
    getColorEvent
from LMT.USV.usvDescription.Compute_Number_USVs import *
'''
import string
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
'''
from LMT.USV.figure.figUtil import getStarsFromPvalues
from LMT.USV.usvDescription.Compute_Number_USVs import getDataFrameWT, getDataFrameKO
from LMT.USV.figure.vocTraitUsagePerEventContext import plotVocTraitUsagePerEventContextPerSetShort
'''

def cleanUSVTimeLine( timeLine , selectedTimeLine ):
    
    #Clean the USV timeline from USV events that were not selected (i.e., where not all animals were detected)
    #note: the function changes the timeline given as argument.
    
    for usvEvent in timeLine.getEventList():
        found = False
        for selectedEvent in selectedTimeLine.getEventList():
            if selectedEvent.startFrame == usvEvent.startFrame and selectedEvent.endFrame == usvEvent.endFrame:
                found = True
                break
                
        if not found:
            timeLine.eventList.remove( usvEvent )

def computeAcousticVariablesPerContext(experiments, tmin, tmax, jsonExtension, eventListToTest, acousticVariables):
    # initialise dictionaries to store data
    data = {}

    for variable in acousticVariables + ['numberOfVoc']:
        data[variable] = {}

        for event in eventListToTest:
            data[variable][event] = {}

            for strain in strainList:
                data[variable][event][strain] = {}

                for sex in sexList:
                    data[variable][event][strain][sex] = {}

                    for age in ageList:
                        data[variable][event][strain][sex][age] = {}

    # extract data from files
    for exp in experiments:
        print('exp: ', exp.getFullName())

        file = exp.file
        print(file)
        connection = sqlite3.connect(file)
        pool = AnimalPool()
        pool.loadAnimals(connection)
        strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile(connection)

        strainFile = strainAgeSexPairGenoFile[0]
        ageFile = strainAgeSexPairGenoFile[1]
        sexFile = strainAgeSexPairGenoFile[2]
        pairFile = strainAgeSexPairGenoFile[3]
        genoFile = strainAgeSexPairGenoFile[4]

        # initialise the datalist for the specific file:
        for variable in acousticVariables + ['numberOfVoc']:
            for event in eventListToTest:
                data[variable][event][strainFile][sexFile][ageFile][pairFile] = []

        # pre load all USVs where both animals are identified.
        print("loading complete usv timeline dictionary")
        USVTimeLine = EventTimeLine(connection, "USV seq", idA=None, minFrame=tmin, maxFrame=tmax,
                                    loadEventIndependently=True)
        selectedUSVTimeLine = EventTimeLine(connection, "Selected USV", idA=None, minFrame=tmin, maxFrame=tmax)

        cleanUSVTimeLine(USVTimeLine, selectedUSVTimeLine)

        print("USV timeline cleaned from unselected USV.")

        # Remove all USV seq where there are some repeated USVs:
        for usvEvent in USVTimeLine.getEventList():
            if usvEvent.metadata['nbBadRepeat'] >= 1:
                USVTimeLine.eventList.remove(usvEvent)

        # clean the voc timeline using the "excluded" metadata
        print("Cleaning voc with excluded metadata")
        vocTimeLine = EventTimeLine(connection, "Voc", minFrame=tmin, maxFrame=tmax, loadEventIndependently=True)
        
        print("remaining voc events: ", len(vocTimeLine.eventList))

        # Select the voc sequences that are overlapping with the event to test:
        for eventToTest in eventListToTest:
            eventToTestTimeLine = EventTimeLine(connection, eventToTest, minFrame=tmin, maxFrame=tmax)
            eventToTestDictionary = eventToTestTimeLine.getDictionary(minFrame=tmin, maxFrame=tmax)

            for usvEvent in USVTimeLine.getEventList():
                for t in range(usvEvent.startFrame, usvEvent.endFrame + 1):
                    if t in eventToTestDictionary:
                        data["numberOfVoc"][eventToTest][strainFile][sexFile][ageFile][pairFile].append(
                            usvEvent.metadata['nbVoc'])
                        break

                        # Select the vocalisations that are overlapping with the event to test:
        for eventToTest in eventListToTest:
            eventToTestTimeLine = EventTimeLine(connection, eventToTest, minFrame=tmin, maxFrame=tmax)
            eventToTestDictionary = eventToTestTimeLine.getDictionary(minFrame=tmin, maxFrame=tmax)
            print("Nb of events {}: {}".format(eventToTest, len(eventToTestTimeLine.getEventList())))

            for vocEvent in vocTimeLine.getEventList():
                for t in range(vocEvent.startFrame, vocEvent.endFrame + 1):
                    if t in eventToTestDictionary:
                        for variable in acousticVariables:
                            data[variable][eventToTest][strainFile][sexFile][ageFile][pairFile].append(
                                vocEvent.metadata[variable])

                        break
        connection.close()

    for variable in acousticVariables:
        for eventToTest in eventListToTest:
            print(eventToTest)
            print("length of dataList: ", len(data[variable][eventToTest][strainFile][sexFile][ageFile][pairFile]))
            print("--------------------")

    # Create a json file to store the computation
    with open("dataAcousticAnalysis" + jsonExtension + ".json", 'w') as fp:
        json.dump(data, fp, indent=4)

    print("json file with acoustic measurements created for ", jsonExtension)
    print("Job done.")


def createDataframeFromJsonVarPerContext( jsonFile, eventListToTest, acousticVariables ):
    #open the json file to work on pre-computed data
    with open( jsonFile ) as json_data:
        dataUsv = json.load(json_data)
    print("json file for acoustic variables re-imported.")
    print(list(dataUsv.keys()))

    dataDic = { 'strain': [], 'sex': [], 'age': [], 'pair': [], 'event': [] }
    for var in acousticVariables:
        dataDic[var] = []

    for event in eventListToTest:
        for strain in strainList:
            for sex in sexList:
                for age in ageList:
                    for pair in list(dataUsv['durationMs'][event][strain][sex][age].keys()):
                        dataList = dataUsv['durationMs'][event][strain][sex][age][pair]
                        print(dataList)
                        print(len(dataList))
                        dataDic['event'].extend( [event] * len(dataList))
                        dataDic['strain'].extend( [strain] * len(dataList))
                        dataDic['sex'].extend ( [sex] * len(dataList) )
                        dataDic['age'].extend( [age] * len(dataList) )
                        dataDic['pair'].extend( [pair] * len(dataList) )
                        for var in acousticVariables:
                            dataDic[var].extend( dataUsv[var][event][strain][sex][age][pair] )

    print('sex: ', len(dataDic['sex']))
    print('pair: ', len(dataDic['pair']))
    print('{}: {}'.format( acousticVariables[0], len(dataDic[acousticVariables[0]])) )
    print( var, len(dataDic[var]))
    df = pd.DataFrame.from_dict(dataDic)
    print('################')
    print(df.head())
    return df


def plotAcousticVarWTBoxplotContext(ax, dataframe, ageClass, eventList, acousticVariables, yMinUsv, yMaxUsv, variable, letter):
    selectedDataframe = dataframe[(dataframe['age'] == ageClass) & (dataframe['strain'] == 'C57BL/6J')]
    yLabel = getFigureLabelTrait(variable)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = list(range(0, len(eventList)))
    ax.set_xticks(xIndex)
    xLabels = []
    for i in eventList:
        xLabels.append(getFigureBehaviouralEventsLabels(i))

    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(eventList)+0.5)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-0.9, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='event', y=variable, hue='sex', data=selectedDataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    bp.legend().set_visible(False)
    if variable == 'durationMs':
        bp.legend(bbox_to_anchor=(1.01, 0.8), borderaxespad=0)
        M_patch = mpatches.Patch(facecolor=getColorAge(ageClass), edgecolor='grey', label='B6 male {}'.format(ageClass))
        F_patch = mpatches.Patch(facecolor=getColorAge(ageClass), edgecolor='black', label='B6 female {}'.format(ageClass))
        bp.legend(handles=[M_patch, F_patch], loc=(0.05, 0.7))
        #bp.legend().set_visible(True)

    colorList = [getColorAge(age=ageClass), getColorAge(age=ageClass)] * len(eventList)
    edgeList = ['grey', 'black'] * len(eventList)
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.set_xticklabels(xLabels, rotation=45, fontsize=12, horizontalalignment='right')

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)
    #bp.legend().set_visible(True)

    k = 0
    for event in eventList:
        valM = selectedDataframe[variable][(selectedDataframe['event'] == event) & (selectedDataframe['sex'] == 'male')]
        valF = selectedDataframe[variable][(selectedDataframe['event'] == event) & (selectedDataframe['sex'] == 'female')]
        nMales = len(valM)
        nFemales = len(valF)
        meanM = np.mean(valM)
        meanF = np.mean(valF)
        print('Tests between sexes: females {} voc, males {} voc'.format(nFemales, nMales))

        ax.text(k - 0.1,
                yMinVar[variable] - 0.02 * (yMaxVar[variable] - yMinVar[variable]),
                'M', rotation=0, fontsize=12, verticalalignment='top', horizontalalignment='center',
                color='black')
        ax.text(k + 0.15,
                yMinVar[variable] - 0.02 * (yMaxVar[variable] - yMinVar[variable]),
                'F', rotation=0, fontsize=12, verticalalignment='top', horizontalalignment='center',
                color='black')

        if variable == 'durationMs':
            ax.text(k - 0.1,
                    yMinVar[variable] - 0.08 * (yMaxVar[variable] - yMinVar[variable]),
                    nMales, rotation=45, fontsize=8, verticalalignment='top', horizontalalignment='right',
                    color='black')
            ax.text(k + 0.15,
                    yMinVar[variable] - 0.08 * (yMaxVar[variable] - yMinVar[variable]),
                    nFemales, rotation=45, fontsize=8, verticalalignment='top', horizontalalignment='right',
                    color='black')

        #tests done only if the number of USVs is large enough (arbitrary threshold) in males
        if nMales >= 20:
            # Mixed model: variable to explain: value; fixed factor = genotype; random effect: pair
            dfTest = dataframe[['age', 'sex', 'pair', 'strain', 'event', variable]]
            dfTest = dfTest[dfTest['event'] == event]
            dfTest.rename(columns={'age': 'age', 'sex': 'sex', 'pair': 'pair', variable: 'value', 'event': 'event'},
                          inplace=True)
            # create model:
            model = smf.mixedlm("value ~ sex", dfTest, groups='pair')
            # run model:
            result = model.fit()
            # print summary
            print('################')
            print(variable)
            print(result.summary())
            p, sign = extractPValueFromLMMResult(result=result, keyword='male')
            '''T, p = ttest_ind( valM, valF )
            print(variable, ' ', ageClass, ' males ', meanM, ' females ', meanF, ' T = ', T, 'p = ', p,
                  getStarsFromPvalues(p, len(eventListToTest) * len(acousticVariables)))'''
            correction = len(eventListToTest) * len(acousticVariables)
            if p * correction >= 0.05:
                stars = getStarsFromPvalues(p, 1)
            elif p * correction < 0.05:
                stars = getStarsFromPvalues(p, 1) + '°'
            ax.text(k,
                    yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]), stars,
                    fontsize=16, horizontalalignment='center', color='black', weight='bold')

        k += 1


def plotAcousticVarKOBoxplotContext(ax, dataframe, eventList, acousticVariables, yMinUsv, yMaxUsv, variable, letter):
    yLabel = getFigureLabelTrait(variable)
    xLabels = []
    for event in eventList:
        xLabels.append(getFigureBehaviouralEventsLabels(event))
    ax.set_xlabel(xLabels)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = list(range(0, len(acousticVariables)))
    ax.set_xticks(xIndex)

    ax.set_ylabel(yLabel, fontsize=18)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(acousticVariables)+1)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='event', y=variable, hue='strain', data=dataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    bp.legend().set_visible(False)


    colorList = [getColorWT(), getColorKO()] * len(eventList)
    edgeList = ['black', 'black']  * len(eventList)
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)
    ax.set_xticklabels(xLabels, rotation=45, fontsize=12, horizontalalignment='center')

    #bp.legend().set_visible(True)

    k = 0
    for event in eventList:
        valWT = dataframe[variable][(dataframe['event'] == event) & (dataframe['strain'] == 'C57BL/6J')]
        valKO = dataframe[variable][(dataframe['event'] == event) & (dataframe['strain'] == 'Shank3')]
        nWT = len(valWT)
        nKO = len(valKO)
        meanWT = np.mean(valWT)
        meanKO = np.mean(valKO)
        print('Tests between sexes: C57BL/6J {} voc, Shank3 {} voc'.format(nWT, nKO))
        if meanWT > meanKO:
            colStar = 'black'
        else:
            colStar = 'red'

        if variable == 'durationMs':
            ax.text(k - 0.1,
                    yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    nWT, rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right',
                    color='black')
            ax.text(k + 0.15,
                    yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    nKO, rotation=45, fontsize=12, verticalalignment='top', horizontalalignment='right',
                    color='black')

        #tests done only if the number of USVs is large enough (arbitrary threshold) in males
        if nWT >= 20:
            # Mixed model: variable to explain: value; fixed factor = genotype; random effect: pair
            dfTest = dataframe[['age', 'sex', 'pair', 'strain', 'event', variable]]
            dfTest = dfTest[dfTest['event']== event]
            dfTest.rename(columns={'age': 'age', 'sex': 'sex', 'pair': 'pair', variable: 'value', 'event': 'event'}, inplace=True)
            # create model:
            model = smf.mixedlm("value ~ strain", dfTest, groups='pair')
            # run model:
            result = model.fit()
            # print summary
            print('################')
            print(variable)
            print(result.summary())
            p, sign = extractPValueFromLMMResult(result=result, keyword='Shank3')
            '''U, p = ttest_ind( valWT, valKO )
            print(variable, ' B6 ', meanWT, ' Shank3 ', meanKO, ' U = ', U, 'p = ', p,
                  getStarsFromPvalues(p, len(eventListToTest) * len(acousticVariables)))'''
            correction = len(eventListToTest) * len(acousticVariables)
            if p * correction >= 0.05:
                stars = getStarsFromPvalues(p, 1)
            elif p * correction < 0.05:
                stars = getStarsFromPvalues(p, 1) + '°'
            ax.text(k,
                    yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]), stars,
                    fontsize=16, horizontalalignment='center', color=colStar, weight='bold')

        k += 1


def plotAcousticVarWTBoxplotContextAge(ax, dataframe, sexClass, eventList, acousticVariables, yMinUsv, yMaxUsv, variable, letter):
    selectedDataframe = dataframe[(dataframe['sex'] == sexClass)]
    yLabel = getFigureLabelTrait(variable)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = list(range(0, len(acousticVariables)))
    ax.set_xticks(xIndex)
    ax.set_xticklabels(acousticVariables, rotation=45, fontsize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=18)
    ax.xaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(acousticVariables)+1)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='event', y=variable, hue='age', data=selectedDataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    bp.legend().set_visible(False)


    colorList = [getColorAge('5we'), getColorAge('3mo'), getColorAge('7mo')] * len(eventList)
    edgeList = ['black', 'black', 'black'] * len(eventList)
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)
    xLabels = []
    for event in eventList:
        xLabels.append(getFigureBehaviouralEventsLabels(event))
    ax.set_xticklabels(xLabels, rotation=45, fontsize=12, horizontalalignment='center')

    #bp.legend().set_visible(True)

    k = 0
    for event in eventList:
        val5we = selectedDataframe[variable][(selectedDataframe['event'] == event) & (selectedDataframe['age'] == '5we')]
        val3mo = selectedDataframe[variable][(selectedDataframe['event'] == event) & (selectedDataframe['age'] == '3mo')]
        val7mo = selectedDataframe[variable][(selectedDataframe['event'] == event) & (selectedDataframe['age'] == '7mo')]
        n5we = len(val5we)
        n3mo = len(val3mo)
        n7mo = len(val7mo)
        mean5we = np.mean(val5we)
        mean3mo = np.mean(val3mo)
        mean7mo = np.mean(val7mo)
        print('Tests between ages: 5we {} voc, 3mo {} voc, 7mo {} voc'.format(n5we, n3mo, n7mo))

        if variable == 'durationMs':
            ax.text(k - 0.15,
                    yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    n5we, rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                    color='black')
            ax.text(k + 0.05,
                    yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    n3mo, rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                    color='black')
            ax.text(k + 0.25,
                    yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    n7mo, rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                    color='black')

        #tests done only if the number of USVs is large enough (arbitrary threshold) in 5we
        if (n5we >= 20) & (n3mo >=20) & (n7mo >= 20):
            U, p = ttest_ind( val5we, val3mo )
            correction = 2 * len(eventListToTest) * len(acousticVariables)
            if p * correction >= 0.05:
                stars = getStarsFromPvalues(p, 1)
            elif p * correction < 0.05:
                stars = getStarsFromPvalues(p, 1) + '°'
            print(variable, ' 5we ', mean5we, ' 3mo ', mean3mo, ' U = ', U, 'p = ', p,
                  stars)
            ax.text(k-0.15, yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]), stars,
                    fontsize=16, horizontalalignment='center', color='black', weight='bold')

            U, p = ttest_ind(val3mo, val7mo)
            correction = 2 * len(eventListToTest) * len(acousticVariables)
            if p * correction >= 0.05:
                stars = getStarsFromPvalues(p, 1)
            elif p * correction < 0.05:
                stars = getStarsFromPvalues(p, 1) + '°'
            print(variable, ' 3mo ', mean3mo, ' 7mo ', mean7mo, ' U = ', U, 'p = ', p,
                  stars)
            ax.text(k+0.15, yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                    stars, fontsize=16, horizontalalignment='center', color='black', weight='bold')

        k += 1


def plotAcousticVarWTBoxplotContextAgeShort(ax, dataframe, sexClass, ageClass, eventList, eventListToTest, acousticVariables, yMinUsv, yMaxUsv, variable, letter):
    if ageClass == None:
        selectedDataframe = dataframe[(dataframe['sex'] == sexClass)]
    else:
        selectedDataframe = dataframe[(dataframe['sex'] == sexClass) & (dataframe['age'] == ageClass)]

    yLabel = variable
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #xIndex = list(range(0, len(acousticVariables)))
    #ax.set_xticks(xIndex)
    #ax.set_xticklabels(eventList, rotation=0, FontSize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, FontSize=18)
    ax.yaxis.set_tick_params(direction="in")
    #ax.set_xlim(0, len(acousticVariables)+1)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='age', y=variable, hue='event', data=selectedDataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)

    y_patch = {}
    patchList = []
    colorListSimple = []

    for n in range(0, len(eventList)):
        y_patch[n] = mpatches.Patch(facecolor=getColorEvent(eventList[n]), edgecolor='black', label=eventList[n], alpha=0.7)
        patchList.append(y_patch[n])
        colorListSimple.append(getColorEvent(eventList[n]))

    bp.legend(loc='center left', bbox_to_anchor=(0.7, 0.5), borderaxespad=0, handles=patchList)

    colorList = colorListSimple * len(ageList)
    edgeList = ['black'] * len(eventList) * len(ageList)
    l = 0
    for box in bp.artists:
        box.set_facecolor(colorList[l])
        box.set_edgecolor(edgeList[l])
        l += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, FontSize=15)
    #bp.legend().set_visible(True)
    k = 0
    if ageClass == None:

        for age in ageList:
            val = {}
            nEvent = {}
            meanEvent = {}
            for n in range(0, len(eventList)):
                val[n] = selectedDataframe[variable][(selectedDataframe['event'] == eventList[n]) & (selectedDataframe['age'] == age)]
                nEvent[n] = len(val[n])
                meanEvent[n] = np.mean(val[n])

            for n in range(0, len(eventList)):
                for m in range(0, len(eventList)):
                    if n < m:
                        print(nEvent[n], meanEvent[n])
                        print(nEvent[m], meanEvent[m])
                        print('k: ', k)
                        print(
                            'Tests between events at {}: {} {} voc, {} {} voc'.format(ageClass, eventList[n], nEvent[n],
                                                                                      eventList[m], nEvent[m]))
                        # tests done only if the number of USVs is large enough (arbitrary threshold)
                        if (nEvent[n] >= 20) & (nEvent[m] >= 20):
                            U, p = ttest_ind(val[n], val[m])
                            print(variable, eventList[n], meanEvent[n], eventList[m], meanEvent[m], ' U = ', U, 'p = ',
                                  p,
                                  getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)))
                            ax.text(k + 0.1 * (n + (m - 1)),
                                    yMaxVar[variable] - 0.06 * (m - n) * (yMaxVar[variable] - yMinVar[variable]),
                                    getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)),
                                    FontSize=16, horizontalalignment='center', color='black', weight='bold')

                if variable == eventList[0]:
                    ax.text(k + 0.15 * (n - 1),
                            yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                            nEvent[n], rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                            color='black')
            k += 1


    else:
        k = 0
        val = {}
        nEvent = {}
        meanEvent = {}
        for n in range(0, len(eventList)):
            val[n] = selectedDataframe[variable][(selectedDataframe['event'] == eventList[n])]
            nEvent[n] = len(val[n])
            meanEvent[n] = np.mean(val[n])

        for n in range(0, len(eventList)):
            for m in range(0, len(eventList)):
                if n < m:
                    print(nEvent[n], meanEvent[n])
                    print(nEvent[m], meanEvent[m])

                    print('Tests between events at {}: {} {} voc, {} {} voc'.format(ageClass, eventList[n], nEvent[n], eventList[m], nEvent[m]))
                    # tests done only if the number of USVs is large enough (arbitrary threshold)
                    if (nEvent[n] >= 20) & (nEvent[m] >= 20):
                        U, p = ttest_ind(val[n], val[m])
                        print(variable, eventList[n], meanEvent[n], eventList[m], meanEvent[m], ' U = ', U, 'p = ', p,
                              getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)))
                        ax.text(k + 0.1* (n+ (m-1)), yMaxVar[variable] - 0.06*(m-n) * (yMaxVar[variable] - yMinVar[variable]),
                                getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)),
                                FontSize=16, horizontalalignment='center', color='black', weight='bold')

            if variable == eventList[0]:
                ax.text(k + 0.15*(n-1),
                        yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                        nEvent[n], rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                        color='black')


def plotAcousticVarKOBoxplotContextShort(ax, dataframe, sexClass, strainClass, eventList, eventListToTest, acousticVariables, yMinUsv, yMaxUsv, variable, letter):

    selectedDataframe = dataframe[(dataframe['sex'] == sexClass)]

    yLabel = variable
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #xIndex = list(range(0, len(acousticVariables)))
    #ax.set_xticks(xIndex)
    #ax.set_xticklabels(eventList, rotation=0, FontSize=12, horizontalalignment='center')
    ax.set_ylabel(yLabel, FontSize=18)
    ax.yaxis.set_tick_params(direction="in")
    #ax.set_xlim(0, len(acousticVariables)+1)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            FontSize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='strain', y=variable, hue='event', data=selectedDataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)

    y_patch = {}
    patchList = []
    colorListSimple = []

    for n in range(0, len(eventList)):
        y_patch[n] = mpatches.Patch(facecolor=getColorEvent(eventList[n]), edgecolor='black', label=eventList[n], alpha=0.7)
        patchList.append(y_patch[n])
        colorListSimple.append(getColorEvent(eventList[n]))

    bp.legend(loc='center left', bbox_to_anchor=(0.7, 0.5), borderaxespad=0, handles=patchList)

    colorList = colorListSimple * len(ageList)
    edgeList = ['black'] * len(eventList) * len(ageList)
    l = 0
    for box in bp.artists:
        box.set_facecolor(colorList[l])
        box.set_edgecolor(edgeList[l])
        l += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, FontSize=15)
    #bp.legend().set_visible(True)

    k = 0
    if strainClass == None:

        for strain in strainList:
            val = {}
            nEvent = {}
            meanEvent = {}
            for n in range(0, len(eventList)):
                val[n] = selectedDataframe[variable][
                    (selectedDataframe['event'] == eventList[n]) & (selectedDataframe['strain'] == strain)]
                nEvent[n] = len(val[n])
                meanEvent[n] = np.mean(val[n])

            for n in range(0, len(eventList)):
                for m in range(0, len(eventList)):
                    if n < m:
                        print(nEvent[n], meanEvent[n])
                        print(nEvent[m], meanEvent[m])
                        print('k: ', k)
                        print(
                            'Tests between events for {}: {} {} voc, {} {} voc'.format(strain, eventList[n], nEvent[n],
                                                                                      eventList[m], nEvent[m]))
                        # tests done only if the number of USVs is large enough (arbitrary threshold)
                        if (nEvent[n] >= 20) & (nEvent[m] >= 20):
                            U, p = ttest_ind(val[n], val[m])
                            print(variable, eventList[n], meanEvent[n], eventList[m], meanEvent[m], ' U = ', U, 'p = ',
                                  p,
                                  getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)))
                            ax.text(k + 0.1 * (n + (m - 1)),
                                    yMaxVar[variable] - 0.06 * (m - n) * (yMaxVar[variable] - yMinVar[variable]),
                                    getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)),
                                    FontSize=16, horizontalalignment='center', color='black', weight='bold')

                if variable == eventList[0]:
                    ax.text(k + 0.15 * (n - 1),
                            yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                            nEvent[n], rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                            color='black')
            k += 1


    else:
        k = 0
        val = {}
        nEvent = {}
        meanEvent = {}
        for n in range(0, len(eventList)):
            val[n] = selectedDataframe[variable][(selectedDataframe['event'] == eventList[n])]
            nEvent[n] = len(val[n])
            meanEvent[n] = np.mean(val[n])

        for n in range(0, len(eventList)):
            for m in range(0, len(eventList)):
                if n < m:
                    print(nEvent[n], meanEvent[n])
                    print(nEvent[m], meanEvent[m])

                    print('Tests between events for {}: {} {} voc, {} {} voc'.format(strainClass, eventList[n], nEvent[n],
                                                                                    eventList[m], nEvent[m]))
                    # tests done only if the number of USVs is large enough (arbitrary threshold)
                    if (nEvent[n] >= 20) & (nEvent[m] >= 20):
                        U, p = ttest_ind(val[n], val[m])
                        print(variable, eventList[n], meanEvent[n], eventList[m], meanEvent[m], ' U = ', U, 'p = ', p,
                              getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)))
                        ax.text(k + 0.1 * (n + (m - 1)),
                                yMaxVar[variable] - 0.06 * (m - n) * (yMaxVar[variable] - yMinVar[variable]),
                                getStarsFromPvalues(p, (len(eventListToTest) * len(eventListToTest))/2 * len(acousticVariables)),
                                FontSize=16, horizontalalignment='center', color='black', weight='bold')

            if variable == eventList[0]:
                ax.text(k + 0.15 * (n - 1),
                        yMinVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]),
                        nEvent[n], rotation=45, FontSize=12, verticalalignment='top', horizontalalignment='right',
                        color='black')


if __name__ == '__main__':
    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec
    rc('font', **{'family': 'serif', 'serif': ['Arial']})
    letterList = list(string.ascii_uppercase)

    while True:
        question = "Do you want to:"
        question += "\n\t [c]ompute the acoustic variables from USVs in each context?"
        question += "\n\t [p]lot each acoustic trait according to the contexts and separated by sex, for each age in B6 mice?"
        question += "\n\t [pko]lot each acoustic trait according to the contexts and separated by genotype?"
        question += "\n\t [page]lot each acoustic trait according to the contexts and separated by age separated by sex in B6?"
        question += "\n\t [pc]lot each acoustic trait according to selected contexts and separated by age in B6 females?"
        question += "\n\t [pc1]lot each acoustic trait according to selected contexts at 5 we in B6 females?"
        question += "\n\t [pc2]lot each acoustic trait according to selected contexts at 3mo in B6 females?"
        question += "\n\t [pko1]lot each acoustic trait according to selected contexts at 3mo in B6 & Shank3 females?"
        question += "\n\t [pko2]lot each acoustic trait of close contacts at 3mo in B6 & Shank3 females?"
        question += "\n\t [pfigbp]plot fig for acoustic variations across contexts in B6 and Shank3 mice?"
        question += "\n\t [pfig]plot fig 5 for acoustic variations across contexts in B6 and Shank3 mice with heatmaps?"
        question += "\n\t [sfig] plot suppl fig with heatmaps for acoustic variations across contexts in B6 females for each age separately?"
        question += "\n\t [sfig1] plot suppl fig with heatmaps for acoustic variations across contexts in Shank3 females?"

        question += "\n"
        answer = input("Action: ")

        if answer == 'c':
            '''
            This codes extracts the acoustic variables of USV according to the contexts of emission. Data are stored in a json file for further computation.
            '''
            print("Code launched.")

            experiments = getAllExperimentList()
            tmin, tmax = getMinTMaxTInput()
            jsonExtension = "_all_pairs"

            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()

            computeAcousticVariablesPerContext(experiments, tmin, tmax, jsonExtension, eventListToTest, acousticVariables)

            break

        if answer == 'p':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariablesComplete = getFigureVocTraits()

            #Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest, acousticVariables=acousticVariablesComplete)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            for age in ageList:
                k = 0
                acousticVariables = acousticVariablesComplete[0:8]
                fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1, figsize=(12, 3*len(acousticVariables)), sharex=True)

                row = 0
                for var in acousticVariables:
                    plotAcousticVarWTBoxplotContext(ax=axes[row], dataframe=dataframeWT, ageClass=age, eventList = eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[k])
                    row += 1
                    k += 1
                plt.tight_layout()
                plt.show()
                fig.savefig('fig_acoustic_var_per_context_sex_{}_1.pdf'.format(age), dpi=300)
                fig.savefig('fig_acoustic_var_per_context_sex_{}_1.jpg'.format(age), dpi=200)

                acousticVariables = acousticVariablesComplete[8:16]
                fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1,
                                         figsize=(12, 3 * len(acousticVariables)), sharex=True)

                row = 0
                for var in acousticVariables:
                    plotAcousticVarWTBoxplotContext(ax=axes[row], dataframe=dataframeWT, ageClass=age,
                                                    eventList=eventListToTest, acousticVariables=acousticVariables,
                                                    yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                    letter=letterList[k])
                    row += 1
                    k += 1
                plt.tight_layout()
                plt.show()
                fig.savefig('fig_acoustic_var_per_context_sex_{}_2.pdf'.format(age), dpi=300)
                fig.savefig('fig_acoustic_var_per_context_sex_{}_2.jpg'.format(age), dpi=200)

            print('Job done.')
            break

        if answer == 'pko':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariablesComplete = getFigureVocTraits()


            #Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest, acousticVariables=acousticVariablesComplete)

            # data wild-type
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')
            k = 0
            acousticVariables = acousticVariablesComplete[0:8]
            fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1, figsize=(12, 3*len(acousticVariables)), sharex=True)

            row = 0
            for var in acousticVariables:
                plotAcousticVarKOBoxplotContext(ax=axes[row], dataframe=dataframeKO, eventList = eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[k])
                row += 1
                k += 1

            ax = axes[0]
            wt_patch = mpatches.Patch(color=getColorWT(), label='C57BL/6J')
            ko_patch = mpatches.Patch(color=getColorKO(), label='Shank3-/-')
            ax.legend(handles=[wt_patch, ko_patch], loc=(0.05, 0.7))

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_acoustic_var_per_context_ko_1.pdf', dpi=300)
            fig.savefig('fig_acoustic_var_per_context_ko_1.jpg', dpi=200)

            acousticVariables = acousticVariablesComplete[8:16]
            fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1, figsize=(12, 3 * len(acousticVariables)),
                                     sharex=True)

            row = 0
            for var in acousticVariables:
                plotAcousticVarKOBoxplotContext(ax=axes[row], dataframe=dataframeKO, eventList=eventListToTest,
                                                acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                                                variable=var, letter=letterList[k])
                row += 1
                k += 1

            ax = axes[0]
            wt_patch = mpatches.Patch(color=getColorWT(), label='C57BL/6J')
            ko_patch = mpatches.Patch(color=getColorKO(), label='Shank3-/-')
            ax.legend(handles=[wt_patch, ko_patch], loc=(0.05, 0.7))

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_acoustic_var_per_context_ko_2.pdf', dpi=300)
            fig.savefig('fig_acoustic_var_per_context_ko_2.jpg', dpi=200)

            print('Job done.')
            break


        if answer == 'page':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariablesComplete = getFigureVocTraits()

            #Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest, acousticVariables=acousticVariablesComplete)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            for sex in sexList:
                acousticVariables = acousticVariablesComplete[0:8]
                k = 0
                fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1, figsize=(12, 3*len(acousticVariables)), sharex=True)

                row = 0
                for var in acousticVariables:
                    plotAcousticVarWTBoxplotContextAge(ax=axes[row], dataframe=dataframeWT, sexClass=sex, eventList = eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[k])
                    row += 1
                    k += 1

                ax = axes[0]
                y_patch = mpatches.Patch(facecolor=getColorAge('5we'), edgecolor='black', label='5 we')
                m_patch = mpatches.Patch(facecolor=getColorAge('3mo'), edgecolor='black', label='3 mo')
                o_patch = mpatches.Patch(facecolor=getColorAge('7mo'), edgecolor='black', label='7 mo')
                ax.legend(handles=[y_patch, m_patch, o_patch], loc=(0.05, 0.6))

                plt.tight_layout()
                plt.show()
                fig.savefig('fig_acoustic_var_per_context_age_{}_1.pdf'.format(sex), dpi=300)
                fig.savefig('fig_acoustic_var_per_context_age_{}_1.jpg'.format(sex), dpi=200)

                acousticVariables = acousticVariablesComplete[8:16]
                fig, axes = plt.subplots(nrows=len(acousticVariables), ncols=1,
                                         figsize=(12, 3 * len(acousticVariables)), sharex=True)

                row = 0
                for var in acousticVariables:
                    plotAcousticVarWTBoxplotContextAge(ax=axes[row], dataframe=dataframeWT, sexClass=sex,
                                                       eventList=eventListToTest, acousticVariables=acousticVariables,
                                                       yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                       letter=letterList[k])
                    row += 1
                    k += 1

                ax = axes[0]
                y_patch = mpatches.Patch(facecolor=getColorAge('5we'), edgecolor='black', label='5 we')
                m_patch = mpatches.Patch(facecolor=getColorAge('3mo'), edgecolor='black', label='3 mo')
                o_patch = mpatches.Patch(facecolor=getColorAge('7mo'), edgecolor='black', label='7 mo')
                ax.legend(handles=[y_patch, m_patch, o_patch], loc=(0.05, 0.7))

                plt.tight_layout()
                plt.show()
                fig.savefig('fig_acoustic_var_per_context_age_{}_2.pdf'.format(sex), dpi=300)
                fig.savefig('fig_acoustic_var_per_context_age_{}_2.jpg'.format(sex), dpi=200)
            print('Job done.')
            break

        if answer == 'pc':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()
            eventListToTestShort = ['Stop isolated', 'Oral-oral Contact', 'Train2']
            acousticVariablesShort = ['durationMs', 'frequencyDynamicHz', 'griffIndex']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            sex = 'female'

            fig, axes = plt.subplots(nrows=1, ncols=len(acousticVariablesShort), figsize=(4*len(acousticVariablesShort), 4), sharex=True)

            row = 0
            for var in acousticVariablesShort:
                plotAcousticVarWTBoxplotContextAgeShort(ax=axes[row], dataframe=dataframeWT, sexClass=sex, ageClass=None, eventList = eventListToTestShort, eventListToTest=eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[row])
                row += 1
            plt.tight_layout()
            plt.show()
            fig.savefig('fig_contexts_{}.pdf'.format(sex), dpi=300)
            print('Job done.')

            break

        if answer == 'pc1':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()
            eventListToTestShort = ['Stop isolated', 'Break contact', 'Oral-oral Contact', 'FollowZone Isolated']
            acousticVariablesShort = ['durationMs', 'frequencyDynamicHz', 'griffIndex']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            sex = 'female'

            fig, axes = plt.subplots(nrows=1, ncols=len(acousticVariablesShort), figsize=(5*len(acousticVariablesShort), 4), sharex=True)

            row = 0
            for var in acousticVariablesShort:
                plotAcousticVarWTBoxplotContextAgeShort(ax=axes[row], dataframe=dataframeWT, sexClass=sex, ageClass='5we', eventList = eventListToTestShort, eventListToTest=eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[row])
                row += 1
            plt.tight_layout()
            plt.show()
            fig.savefig('fig_contexts_{}_5we.pdf'.format(sex), dpi=300)
            print('Job done.')

            break

        if answer == 'pc2':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()
            eventListToTestShort = ['Oral-oral Contact', 'Oral-genital Contact', 'Side by side Contact']
            acousticVariablesShort = ['frequencyDynamicHz']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            sex = 'female'

            fig, axes = plt.subplots(nrows=1, ncols=len(acousticVariablesShort), figsize=(5*len(acousticVariablesShort), 4), sharex=True)

            row = 0
            for var in acousticVariablesShort:
                if len(acousticVariablesShort) == 1:
                    ax = axes
                else:
                    ax = axes[row]
                plotAcousticVarWTBoxplotContextAgeShort(ax=ax, dataframe=dataframeWT, sexClass=sex, ageClass='3mo', eventList = eventListToTestShort, eventListToTest=eventListToTest, acousticVariables=acousticVariables, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[row])
                row += 1
            plt.tight_layout()
            plt.show()
            fig.savefig('fig_contexts_{}_3mo.pdf'.format(sex), dpi=300)
            print('Job done.')

            break


        if answer == 'pko1':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()
            eventListToTestShort = ['Oral-genital Contact', 'FollowZone Isolated', 'Train2']
            acousticVariablesShort = ['durationMs', 'frequencyDynamicHz', 'griffIndex']
            sex = 'female'
            #Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort, acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')

            fig, axes = plt.subplots(nrows=1, ncols=len(acousticVariablesShort), figsize=(5 * len(acousticVariablesShort), 4), sharex=True)

            row = 0
            for var in acousticVariablesShort:
                if len(acousticVariablesShort) == 1:
                    ax = axes
                else:
                    ax = axes[row]
                plotAcousticVarKOBoxplotContextShort(ax=ax, dataframe=dataframeKO, sexClass=sex, strainClass=None,
                                                        eventList=eventListToTestShort, eventListToTest=eventListToTest,
                                                        acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                        yMaxUsv=yMaxVar, variable=var, letter=letterList[row])
                row += 1
            plt.tight_layout()
            plt.show()
            fig.savefig('fig_contexts_{}_ko.pdf'.format(sex), dpi=300)
            print('Job done.')

            break

        if answer == 'pko2':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()
            eventListToTestShort = ['Oral-oral Contact', 'Oral-genital Contact', 'Side by side Contact', 'Side by side Contact, opposite way']
            acousticVariablesShort = ['durationMs', 'frequencyDynamicHz', 'griffIndex']
            sex = 'female'
            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')
            selectedDataframeKO = dataframeKO[(dataframeKO['strain'] == 'Shank3')]

            fig, axes = plt.subplots(nrows=1, ncols=len(acousticVariablesShort),
                                     figsize=(5 * len(acousticVariablesShort), 4), sharex=True)

            row = 0
            for var in acousticVariablesShort:
                if len(acousticVariablesShort) == 1:
                    ax = axes
                else:
                    ax = axes[row]
                plotAcousticVarKOBoxplotContextShort(ax=ax, dataframe=selectedDataframeKO, sexClass=sex, strainClass='Shank3',
                                                     eventList=eventListToTestShort, eventListToTest=eventListToTest,
                                                     acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                     yMaxUsv=yMaxVar, variable=var, letter=letterList[row])
                row += 1
            plt.tight_layout()
            plt.show()
            fig.savefig('fig_contexts_{}_ko_close_contacts.pdf'.format(sex), dpi=300)
            print('Job done.')

            break

        if answer == 'pfigbp':
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()

            # Create general figure
            gs = gridspec.GridSpec(3, 3)
            fig = plt.figure(figsize=(12, 9))
            letter = 0

            # acoustic variations at 3 months between sexes ########################################################
            ax = fig.add_subplot(gs[0, 0:3])
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = ['durationMs']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest,
                                                      acousticVariables=acousticVariables)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            age = '3mo'

            row = 0
            for var in acousticVariables:
                plotAcousticVarWTBoxplotContext(ax=ax, dataframe=dataframeWT, ageClass=age,
                                                eventList=eventListToTest, acousticVariables=acousticVariables,
                                                yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                letter=letterList[letter])
                row += 1

            letter += 1

            #acoustic variations at 5 weeks ########################################################
            ax = fig.add_subplot(gs[1, 0:1])

            # Create dataframe from json file:
            eventListToTestShort = ['Stop isolated', 'Break contact', 'Oral-oral Contact', 'FollowZone Isolated']
            acousticVariablesShort = ['durationMs']
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            sex = 'female'
            for var in acousticVariablesShort:
                plotAcousticVarWTBoxplotContextAgeShort(ax=ax, dataframe=dataframeWT, sexClass=sex,
                                                            ageClass='5we', eventList=eventListToTestShort,
                                                            eventListToTest=eventListToTest,
                                                            acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                            yMaxUsv=yMaxVar, variable=var, letter=letterList[letter])

            letter += 1

            # acoustic variations at 3 months ########################################################
            ax = fig.add_subplot(gs[1, 1:2])
            eventListToTestShort = ['Oral-oral Contact', 'Oral-genital Contact', 'Side by side Contact']
            acousticVariablesShort = ['frequencyDynamicHz']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            sex = 'female'
            for var in acousticVariablesShort:
                plotAcousticVarWTBoxplotContextAgeShort(ax=ax, dataframe=dataframeWT, sexClass=sex, ageClass='3mo',
                                                            eventList=eventListToTestShort, eventListToTest=eventListToTest,
                                                            acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                            yMaxUsv=yMaxVar, variable=var, letter=letterList[letter])

            letter += 1

            # acoustic variations in shank3 follow ########################################################
            ax = fig.add_subplot(gs[1, 2:3])
            eventListToTestShort = ['Oral-genital Contact', 'FollowZone Isolated', 'Train2']
            acousticVariablesShort = ['durationMs']
            sex = 'female'
            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data KO
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')
            for var in acousticVariablesShort:
                plotAcousticVarKOBoxplotContextShort(ax=ax, dataframe=dataframeKO, sexClass=sex, strainClass=None,
                                                         eventList=eventListToTestShort, eventListToTest=eventListToTest,
                                                         acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                         yMaxUsv=yMaxVar, variable=var, letter=letterList[letter])
            letter += 1

            # acoustic variations in shank3 close contacts ########################################################
            eventListToTestShort = ['Oral-oral Contact', 'Oral-genital Contact', 'Side by side Contact',
                                    'Side by side Contact, opposite way']
            acousticVariablesShort = ['durationMs', 'frequencyDynamicHz', 'griffIndex']
            sex = 'female'
            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTestShort,
                                                      acousticVariables=acousticVariablesShort)

            # data Shank3
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')
            selectedDataframeKO = dataframeKO[(dataframeKO['strain'] == 'Shank3')]

            col = 0
            for var in acousticVariablesShort:
                ax = fig.add_subplot(gs[2, col])
                plotAcousticVarKOBoxplotContextShort(ax=ax, dataframe=selectedDataframeKO, sexClass=sex,
                                                     strainClass='Shank3',
                                                     eventList=eventListToTestShort, eventListToTest=eventListToTest,
                                                     acousticVariables=acousticVariables, yMinUsv=yMinVar,
                                                     yMaxUsv=yMaxVar, variable=var, letter=letterList[letter])
                letter += 1
                col += 1




            plt.tight_layout()
            plt.show()
            fig.savefig('fig_acoustic_variations_contexts.pdf', dpi=300)
            print('Job done.')

            break

        if answer == "pfiglegend":
            featureHeatMapPValLegend()
            featureHeatMapEffectSizeLegend()

        if answer == 'pfig':
            '''This script provides Figure 5'''
            print('Code launched.')
            jsonFile = 'dataAcousticAnalysis_all_pairs.json'
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            acousticVariables = getFigureVocTraits()

            # Create general figure
            gs = gridspec.GridSpec(3, 3)
            fig = plt.figure(figsize=(12, 12))
            letter = 0

            '''# acoustic variations at 3 months between sexes ########################################################
            ax = fig.add_subplot(gs[0, 0:3])
            acousticVariables = ['durationMs']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest,
                                                      acousticVariables=acousticVariables)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            age = '3mo'

            row = 0
            for var in acousticVariables:
                plotAcousticVarWTBoxplotContext(ax=ax, dataframe=dataframeWT, ageClass=age,
                                                eventList=eventListToTest, acousticVariables=acousticVariables,
                                                yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                letter=letterList[letter])
                row += 1

            letter += 1
'''
            # acoustic variations at 3 months between sexes in one specific context ########################################################
            ax = fig.add_subplot(gs[0, 0:1])
            acousticVariables = ['durationMs']
            eventListToTest = ['Stop isolated']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest,
                                                      acousticVariables=acousticVariables)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            age = '3mo'

            for var in acousticVariables:
                plotAcousticVarWTBoxplotContext(ax=ax, dataframe=dataframeWT, ageClass=age,
                                                eventList=eventListToTest, acousticVariables=acousticVariables,
                                                yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                letter=letterList[letter])

            letter += 1

            # acoustic variations at 3 months between sexes in one specific context ########################################################
            ax = fig.add_subplot(gs[0, 1:2])
            acousticVariables = ['nbModulation']
            eventListToTest = ['Stop isolated']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest,
                                                      acousticVariables=acousticVariables)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            age = '3mo'

            for var in acousticVariables:
                plotAcousticVarWTBoxplotContext(ax=ax, dataframe=dataframeWT, ageClass=age,
                                                eventList=eventListToTest, acousticVariables=acousticVariables,
                                                yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                letter=letterList[letter])

            letter += 1

            # acoustic variations at 3 months between sexes in one specific context ########################################################
            ax = fig.add_subplot(gs[0, 2:3])
            acousticVariables = ['meanFrequencyHz']
            eventListToTest = ['Stop isolated']

            # Create dataframe from json file:
            df = createDataframeFromJsonVarPerContext(jsonFile=jsonFile, eventListToTest=eventListToTest,
                                                      acousticVariables=acousticVariables)

            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            age = '3mo'

            for var in acousticVariables:
                plotAcousticVarWTBoxplotContext(ax=ax, dataframe=dataframeWT, ageClass=age,
                                                eventList=eventListToTest, acousticVariables=acousticVariables,
                                                yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var,
                                                letter=letterList[letter])

            letter += 1

            #################################################################################
            # acoustic variations between contexts at 3mo in B6 females with heatmaps ########################################################
            experiments = getExperimentList(age="3mo", sex="female", genotype="WT")
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            eventListLabels = []
            for i in eventListToTest:
                eventListLabels.append(getFigureBehaviouralEventsLabels(i))

            acousticVariables = ['durationMs', 'nbModulation', 'meanFrequencyHz']

            tickPos = []
            for i in list(range(len(eventListLabels))):
                tickPos.append(i+0.5)

            #vocTraitUsagePerEventContext.py
            ax = fig.add_subplot(gs[1, 0:1])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            ax.text(-6, 5, 'C57BL/6J female 3mo', fontsize=16, rotation=90, verticalalignment='center', horizontalalignment='center', color=getColorWT(),
                    weight='bold')
            variable = acousticVariables[0]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json', experiments=experiments, genotype='WT', eventListToTest = eventListToTest, variable=variable, ax=ax)
            g.set_yticks(tickPos)
            g.set_yticklabels(eventListLabels, rotation=0, fontsize=10)
            #g.set_xticks(tickPos)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)

            '''# add scale on the heatmaps
            image = 'scale.jpg'
            imgPos = (4.5, -0.5)
            behavSchema = mpimg.imread(image)
            imgBox = OffsetImage(behavSchema, zoom=0.15)
            imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
            g.add_artist(imageBox)
            '''
            letter += 1


            ax = fig.add_subplot(gs[1, 1:2])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            variable = acousticVariables[1]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json', experiments=experiments, genotype='WT',
                                                        eventListToTest=eventListToTest, variable=variable, ax=ax)
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)
            #g.set_xticks(tickPos)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
            letter += 1

            ax = fig.add_subplot(gs[1, 2:3])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            variable = acousticVariables[2]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json', experiments=experiments, genotype='WT',
                                                        eventListToTest=eventListToTest, variable=variable, ax=ax)
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)
            #g.set_xticks(tickPos)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
            letter += 1

            #################################################################################
            # acoustic variations between contexts at 3mo in Shank3 females with heatmaps ########################################################
            experiments = getExperimentList(age="3mo", sex="female", genotype="KO")
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            eventListLabels = []
            for i in eventListToTest:
                eventListLabels.append(getFigureBehaviouralEventsLabels(i))

            #acousticVariables = ['durationMs', 'frequencyDynamicHz', 'griffIndex']

            # vocTraitUsagePerEventContext.py
            ax = fig.add_subplot(gs[2, 0:1])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            ax.text(-6, 5, 'Shank3-/- 3mo', fontsize=16, rotation=90, verticalalignment='center', horizontalalignment='center', color=getColorKO(),
                    weight='bold')

            variable = acousticVariables[0]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json',
                                                            experiments=experiments, genotype='KO',
                                                            eventListToTest=eventListToTest, variable=variable, ax=ax)
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)
            #g.set_xticks(tickPos)
            g.set_yticks(tickPos)
            g.set_yticklabels(eventListLabels, rotation=0, fontsize=10)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')


            letter += 1

            ax = fig.add_subplot(gs[2, 1:2])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            variable = acousticVariables[1]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json',
                                                        experiments=experiments, genotype='KO',
                                                        eventListToTest=eventListToTest, variable=variable, ax=ax)
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)
            #g.set_xticks(tickPos)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
            letter += 1

            ax = fig.add_subplot(gs[2, 2:3])
            ax.text(-1.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black', weight='bold')
            variable = acousticVariables[2]
            g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json',
                                                        experiments=experiments, genotype='KO',
                                                        eventListToTest=eventListToTest, variable=variable, ax=ax)
            ax.set_title(getFigureLabelTrait(variable), fontsize=14)
            #g.set_xticks(tickPos)
            g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
            letter += 1

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_5_acoustic_variations_contexts_heatmaps.pdf', dpi=300)
            fig.savefig('fig_5_acoustic_variations_contexts_heatmaps.jpg', dpi=300)
            print('Job done.')


            break

        if answer == 'sfig':
            #plot the supplementary figures with heatmaps
            acousticVariables = getFigureVocTraits()
            for age in ageList:
                letter = 0
                row = 0
                col = 0
                fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(16, 16))
                # acoustic variations between contexts at 3mo in B6 females with heatmaps ########################################################
                experiments = getExperimentList(age=age, sex="female", genotype="WT")
                eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
                eventListLabels = []
                for i in eventListToTest:
                    eventListLabels.append(getFigureBehaviouralEventsLabels(i))

                tickPos = []
                for i in list(range(len(eventListLabels))):
                    tickPos.append(i + 0.5)

                for variable in acousticVariables:
                    # vocTraitUsagePerEventContext.py
                    ax = axes[row][col]
                    ax.text(-3.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black',
                            weight='bold')

                    g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json',
                                                                     experiments=experiments, genotype='WT',
                                                                     eventListToTest=eventListToTest, variable=variable, ax=ax)
                    g.set_yticks(tickPos)
                    g.set_yticklabels(eventListLabels, rotation=0, fontsize=10)
                    # g.set_xticks(tickPos)
                    g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
                    ax.set_title(getFigureLabelTrait(variable), fontsize=14)
                    if col < 3 :
                        col += 1
                    else:
                        col = 0
                        row += 1
                    letter += 1

                plt.tight_layout()
                plt.show()
                fig.savefig('fig_suppl_heatmaps_acoustic_variations_contexts_{}.pdf'.format(age), dpi=300)
                fig.savefig('fig_suppl_heatmaps_acoustic_variations_contexts_{}.jpg'.format(age), dpi=300)

            print('Job done.')
            break

        if answer == 'sfig1':
            #plot the supplementary figures with heatmaps
            acousticVariables = getFigureVocTraits()
            age = '3mo'
            letter = 0
            row = 0
            col = 0
            fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(16, 16))
            # acoustic variations between contexts at 3mo in B6 females with heatmaps ########################################################
            experiments = getExperimentList(age=age, sex="female", genotype="KO")
            eventListToTest = getFigureBehaviouralEvents(longList=False, withUrinate=False)
            eventListLabels = []
            for i in eventListToTest:
                eventListLabels.append(getFigureBehaviouralEventsLabels(i))

            tickPos = []
            for i in list(range(len(eventListLabels))):
                tickPos.append(i + 0.5)

            for variable in acousticVariables:
                # vocTraitUsagePerEventContext.py
                ax = axes[row][col]
                ax.text(-3.2, -1, letterList[letter], fontsize=20, horizontalalignment='center', color='black',
                        weight='bold')

                g = plotVocTraitUsagePerEventContextPerSetShort2(jsonFile='vocTraitUsagePerEventContext.json',
                                                                 experiments=experiments, genotype='KO',
                                                                 eventListToTest=eventListToTest, variable=variable, ax=ax)
                g.set_yticks(tickPos)
                g.set_yticklabels(eventListLabels, rotation=0, fontsize=10)
                # g.set_xticks(tickPos)
                g.set_xticklabels(eventListLabels, rotation=45, fontsize=10, horizontalalignment='right')
                ax.set_title(getFigureLabelTrait(variable), fontsize=14)
                if col < 3 :
                    col += 1
                else:
                    col = 0
                    row += 1
                letter += 1

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_suppl_heatmaps_acoustic_variations_contexts_Shank3.pdf', dpi=300)
            fig.savefig('fig_suppl_heatmaps_acoustic_variations_contexts_Shank3.jpg', dpi=300)

            print('Job done.')
            break


