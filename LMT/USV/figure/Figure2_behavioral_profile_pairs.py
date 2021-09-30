'''
Created by E. Ey on 07/12/2020
'''

from lmtanalysis.Event import *
from lmtanalysis.Measure import *
import numpy as np;

from scripts.ComputeMeasuresIdentityProfileOneMouseAutomatic import mergeProfileOverNights, getProfileValues, \
    getProfileValuesPairs

np.random.seed(0)
from tkinter.filedialog import askopenfilename
from lmtanalysis.Util import getMinTMaxTAndFileNameInput, getMinTMaxTInput
import sqlite3
from lmtanalysis.FileUtil import getFilesToProcess, getJsonFileToProcess, getStarsFromPvalues
from lmtanalysis.Animal import AnimalPool
from collections import Counter
from LMT.USV.lib.vocUtil import *
import pandas as pd
import seaborn as sns
#import pingouin as pg
from scipy import stats
from scipy.stats.stats import spearmanr, mannwhitneyu
from statsmodels.stats.anova import AnovaRM
from LMT.USV.experimentList.experimentList import getExperimentList,\
    getAllExperimentList
from scipy.stats.morestats import wilcoxon
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from LMT.USV.burster.burster import *
from USV.figure.figUtil import addJitter
from USV.figure.figParameter import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, kruskal, ttest_ind
import string
import matplotlib.patches as mpatches


def singlePlotPerEventProfileBothSexesPerAge(dataframe, valueCat, behavEvent, ax, letter,
                                             xPos, yMin, yMax):
    # image, imgPos, zoom,
    if behavEvent != 'totalDistance':
        event = behavEvent + valueCat

    elif behavEvent == 'totalDistance':
        event = behavEvent
    print("event: ", event)

    ax.text(x=-1, y = yMax + 0.1 * (yMax - yMin), s=letter, fontsize=20, horizontalalignment='center', color='black',
            weight='bold')

    df = dataframe.loc[(dataframe['strain'] == 'C57BL/6J'), :]

    if (valueCat == ' TotalLen'):
        y = [yval/30 for yval in df['value']]
    else:
        y = df['value']

    bp = sns.boxplot(x=df['age'], y=y, hue=df['sex'], hue_order=['male', 'female'], ax=ax, linewidth=1,
                     showmeans=True,
                     meanprops={"marker": 'o',
                                "markerfacecolor": 'white',
                                "markeredgecolor": 'black',
                                "markersize": '8'}, showfliers=False, width=0.8, dodge=True)
    # Change the colors of the border of the boxes
    edgeColor = ['grey', 'black', 'grey', 'black', 'grey', 'black']
    faceColor = [getColorAge('5we'), getColorAge('5we'), getColorAge('3mo'), getColorAge('3mo'), getColorAge('7mo'), getColorAge('7mo')]
    k = 0
    for patch in bp.artists:
        patch.set_edgecolor(edgeColor[k])
        patch.set_facecolor(faceColor[k])
        patch.set_linewidth(1.2)
        k += 1

    #sns.stripplot(x=df['age'], y=df['value'], hue=df['sex'], hue_order=['male', 'female'], jitter=True, color='black', s=5, dodge=True, ax=ax)

    #plot points for each pair across age classes
    groupLevels = list(Counter(df['group']))
    markerList = ['o', 'v', '^', 's', '*', 'h', 'P', 'd']
    m = 0
    for group in groupLevels:
        print(group)
        gdf = df.loc[(df['group'] == group), :]
        if (valueCat == ' TotalLen'):
            y = [yval / 30 for yval in gdf['value']]
        else:
            y = gdf['value']
        sns.stripplot(x=gdf['age'], y=y, hue=gdf['sex'], hue_order=['male', 'female'], jitter=True, color='black',
                  s=5, dodge=True, ax=ax, marker=markerList[m])
        m += 1

    ax.set_ylim(ymin=yMin, ymax=yMax)
    ax.set_xlim(xmin=-0.5, xmax=2.5)
    ax.set_xticks(xPos)
    #ax.set_title(behavEvent, y=1, fontsize=14)
    ax.xaxis.set_tick_params(direction="in")
    ax.tick_params(axis='x', labelsize=14)
    ax.yaxis.set_tick_params(direction="in")
    ax.tick_params(axis='y', labelsize=12)
    if valueCat == ' TotalLen':
        ylabel = 'total duration (s)'
    if valueCat == ' Nb':
        ylabel = 'occurrences'
    if valueCat == ' MeanDur':
        ylabel = 'mean duration (s)'
    elif valueCat == '':
        ylabel = event + ' (m)'
    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def testProfileDataBetweenSexesWT(ageList, df, ax, event, yMinDic, yMaxDic):
    k = 0
    for age in ageList:
        valMale = df['value'][(df['sex'] == 'male') & (df['age'] == age) & (df['strain'] == 'C57BL/6J')]
        valFemale = df['value'][(df['sex'] == 'female') & (df['age'] == age) & (df['strain'] == 'C57BL/6J')]
        print(valMale)
        print(valFemale)
        W, p = mannwhitneyu(valMale, valFemale)
        print('Mann-Whitney U test for {}: W={} p={}'.format(event, W, p))
        correction = len(ageList)
        if p*correction >= 0.05:
            stars = getStarsFromPvalues(p, 1)
        elif p*correction < 0.05:
            stars = getStarsFromPvalues(p, 1)+'°'
        ax.text(x=xPos[k], y=yMinDic[event] + 0.05 * (yMaxDic[event] - yMinDic[event]),
                s=stars, fontsize=16, horizontalalignment='center', color='black',
                weight='bold')
        k += 1

def testProfileDataBetweenAgesWT(ageList, df, ax, event, yMinDic, yMaxDic):
    colorSex = {'male': 'grey', 'female': 'black'}
    xPosSex = {'male': -0.25, 'female': 0.15}
    yPosSex = {'male': 2, 'female': 1}
    for sex in ['male', 'female']:
        adf = df.loc[(df['strain'] == 'C57BL/6J') & (df['sex']==sex), :]
        print(adf)
        ageData = {}
        groupData = {}
        for age in ageList:
            ageData[age] = adf.loc[(adf['age']==age), 'value']
            groupData[age] = adf.loc[(adf['age']==age), 'group']

        print('data 5we: ', ageData['5we'])
        print('group5we: ', groupData['5we'])
        print('group3mo: ', groupData['3mo'])
        print('group3mo: ', groupData['3mo'])
        K, pk = kruskal(ageData['5we'], ageData['3mo'], ageData['7mo'])
        print('pk ', pk)
        if (pk < 0.05):
            W1, p1 = wilcoxon(ageData['5we'], ageData['3mo'])
            print('Wilcoxon U test for {}: W={} p={}'.format(event, W1, p1))
            correction = 2
            if p1 * correction >= 0.05:
                stars = getStarsFromPvalues(p1, 1)
            elif p1 * correction < 0.05:
                stars = getStarsFromPvalues(p1, 1) + '°'
            ax.text(x=0.5+xPosSex[sex], y=yMaxDic[event] - 0.07 * yPosSex[sex] * (yMaxDic[event] - yMinDic[event]),
                    s=stars, fontsize=16, verticalalignment='center', horizontalalignment='center', color=colorSex[sex], weight='bold')
            ax.axhline(y=yMaxDic[event] - 0.08 * yPosSex[sex] * (yMaxDic[event] - yMinDic[event]), xmin=0.16+xPosSex[sex]/6, xmax=0.5+xPosSex[sex]/6, color=colorSex[sex], linewidth=0.5)
            W2, p2 = wilcoxon(ageData['3mo'], ageData['7mo'])
            print('Wilcoxon U test for {}: W={} p={}'.format(event, W2, p2))
            correction = 2 * len(getFigureBehaviouralEvents(longList=False, withUrinate=False))
            if p2 * correction >= 0.05:
                stars = getStarsFromPvalues(p2, 1)
            elif p2 * correction < 0.05:
                stars = getStarsFromPvalues(p2, 1) + '°'
            ax.text(x=1.5+xPosSex[sex], y=yMaxDic[event] - 0.08 * yPosSex[sex] * (yMaxDic[event] - yMinDic[event]),
                    s=stars, fontsize=16, horizontalalignment='center', verticalalignment='center', color=colorSex[sex], weight='bold')
            ax.axhline(y=yMaxDic[event] - 0.11 * yPosSex[sex] * (yMaxDic[event] - yMinDic[event]),
                       xmin=0.5+xPosSex[sex]/6 , xmax=0.85+xPosSex[sex]/6 , color=colorSex[sex], linewidth=0.5)


def singlePlotPerEventProfilePerGenotype(dataframe, valueCat, behavEvent, ax, letter,
                                             xPos, yMin, yMax):
    # image, imgPos, zoom,
    if behavEvent != 'totalDistance':
        event = behavEvent + valueCat

    elif behavEvent == 'totalDistance':
        event = behavEvent
    print("event: ", event)

    ax.text(x=-0.75, y = yMax + 0.1 * (yMax - yMin), s=letter, fontsize=20, horizontalalignment='center', color='black',
            weight='bold')

    df = dataframe.loc[(dataframe['sex'] == 'female') & (dataframe['age'] == '3mo'), :]

    if (valueCat == ' TotalLen'):
        y = [yval/30 for yval in df['value']]
    else:
        y = df['value']

    bp = sns.boxplot(x=df['age'], y=y, hue=df['strain'], hue_order=['C57BL/6J', 'Shank3'], ax=ax, linewidth=1,
                     showmeans=True,
                     meanprops={"marker": 'o',
                                "markerfacecolor": 'white',
                                "markeredgecolor": 'black',
                                "markersize": '8'}, showfliers=False, width=0.4, dodge=True)
    # Change the colors of the border of the boxes
    edgeColor = ['black', 'black']
    faceColor = [getColorWT(), getColorKO()]
    k = 0
    for patch in bp.artists:
        patch.set_edgecolor(edgeColor[k])
        patch.set_facecolor(faceColor[k])
        patch.set_linewidth(1.2)
        k += 1

    #sns.stripplot(x=df['age'], y=df['value'], hue=df['sex'], hue_order=['male', 'female'], jitter=True, color='black', s=5, dodge=True, ax=ax)

    #plot points for each pair across age classes
    groupLevels = list(Counter(df['group']))
    markerList = ['o', 'v', '^', 's', '*', 'h', 'P', 'd', 'o', 'v', '^', 's']
    m = 0
    for group in groupLevels:
        print(group)
        gdf = df.loc[(df['group'] == group), :]
        if (valueCat == ' TotalLen'):
            y = [yval / 30 for yval in gdf['value']]
        else:
            y = gdf['value']
        sns.stripplot(x=gdf['age'], y=y, hue=gdf['strain'], hue_order=['C57BL/6J', 'Shank3'], jitter=True, color='black',
                  s=5, dodge=True, ax=ax, marker=markerList[m])
        m += 1

    ax.set_ylim(ymin=yMin, ymax=yMax)
    ax.set_xlim(xmin=-0.5, xmax=0.5)
    ax.set_xticks(xPos)
    #ax.set_title(behavEvent, y=1, fontsize=14)
    ax.xaxis.set_tick_params(direction="in")
    ax.tick_params(axis='x', labelsize=14)
    #ax.set_xticklabels('3mo')
    ax.yaxis.set_tick_params(direction="in")
    ax.tick_params(axis='y', labelsize=12)
    if valueCat == ' TotalLen':
        ylabel = 'total duration (s)'
    if valueCat == ' Nb':
        ylabel = 'occurrences'
    if valueCat == ' MeanDur':
        ylabel = 'mean duration (s)'
    elif valueCat == '':
        ylabel = event + ' (m)'
    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend().set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

def testProfileDataBetweenGenotypes(df, ax, event, yMinDic, yMaxDic):
    k = 0
    valWt = df['value'][(df['sex'] == 'female') & (df['age'] == '3mo') & (df['strain'] == 'C57BL/6J')]
    valKo = df['value'][(df['sex'] == 'female') & (df['age'] == '3mo') & (df['strain'] == 'Shank3')]
    print(valWt)
    print(valKo)
    W, p = mannwhitneyu(valWt, valKo)
    print('Mann-Whitney U test for {}: W={} p={}'.format(event, W, p))
    correction = 1
    if p * correction >= 0.05:
        stars = getStarsFromPvalues(p, 1)
    elif p * correction < 0.05:
        stars = getStarsFromPvalues(p, 1) + '°'
    ax.text(x=xPos[k], y=yMaxDic[event] - 0.1 * (yMaxDic[event] - yMinDic[event]),
            s=stars, fontsize=16, horizontalalignment='center', color='black',
            weight='bold')



if __name__ == '__main__':

    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})

    behaviouralEventOneMouse = ["Move isolated", "Move in contact", "WallJump", "Stop isolated", "Rear isolated",
                                "Rear in contact", "Contact",
                                #"Group2", "Group3",
                                "Oral-oral Contact",
                                "Oral-genital Contact", "Side by side Contact", "Side by side Contact, opposite way", "Train2", "FollowZone Isolated", "Social approach",
                                "Approach contact", "Get away", "Break contact",
                                #"Group 3 make", "Group 4 make",
                                #"Group 3 break", "Group 4 break",
                                "totalDistance", "experiment"
                                ]

    behaviouralEventOneMouseShort = ["Move isolated", "Stop isolated", "Contact",
                                "Oral-oral Contact", "Oral-genital Contact", "Side by side Contact", "Side by side Contact, opposite way",
                                "Train2", "FollowZone Isolated", "Social approach",
                                "Approach contact", "Get away", "Break contact"
                                ]

    ageList = ['5we', '3mo', '7mo']
    yMinDic = {'totalDistance': 2000,
               "Move isolated TotalLen": 0, "Move in contact TotalLen": 0, "WallJump TotalLen": 0,
               "Stop isolated TotalLen": 30000, "Rear isolated TotalLen": 0,
               "Rear in contact TotalLen": 0, "Contact TotalLen": 000, "Group2 TotalLen": 0,
               "Group3 TotalLen": 0, "Oral-oral Contact TotalLen": 0,
               "Oral-genital Contact TotalLen": 0, "Side by side Contact TotalLen": 0,
               "Side by side Contact, opposite way TotalLen": 0, "Train2 TotalLen": 0,
               "FollowZone Isolated TotalLen": 0, "Social approach TotalLen": 000,
               "Approach contact TotalLen": 00, "Get away TotalLen": 00, "Break contact TotalLen": 00,

               "Move isolated MeanDur": 0, "Move in contact MeanDur": 0, "WallJump MeanDur": 0,
               "Stop isolated MeanDur": 0, "Rear isolated MeanDur": 0,
               "Rear in contact MeanDur": 0, "Contact MeanDur": 0, "Group2 MeanDur": 0,
               "Group3 MeanDur": 0, "Oral-oral Contact MeanDur": 0,
               "Oral-genital Contact MeanDur": 0, "Side by side Contact MeanDur": 0,
               "Side by side Contact, opposite way MeanDur": 0, "Train2 MeanDur": 0,
               "FollowZone Isolated MeanDur": 0, "Social approach MeanDur": 0,
               "Approach contact MeanDur": 0, "Get away MeanDur": 0, "Break contact MeanDur": 0,

               "Move isolated Nb": 0, "Move in contact Nb": 0, "WallJump Nb": 0,
               "Stop isolated Nb": 0, "Rear isolated Nb": 0,
               "Rear in contact Nb": 0, "Contact Nb": 0, "Group2 Nb": 0,
               "Group3 Nb": 0, "Oral-oral Contact Nb": 0,
               "Oral-genital Contact Nb": 0, "Side by side Contact Nb": 0,
               "Side by side Contact, opposite way Nb": 0, "Train2 Nb": 0,
               "FollowZone Isolated Nb": 0, "Social approach Nb": 0,
               "Approach contact Nb": 0, "Get away Nb": 0, "Break contact Nb": 0
               }

    yMaxDic = {'totalDistance': 7000,
               "Move isolated TotalLen": 20000, "Move in contact TotalLen": 20000, "WallJump TotalLen": 400,
               "Stop isolated TotalLen": 100000, "Rear isolated TotalLen": 7000,
               "Rear in contact TotalLen": 4000, "Contact TotalLen": 40000, "Group2 TotalLen": 40000,
               "Group3 TotalLen": 40000, "Oral-oral Contact TotalLen": 7000,
               "Oral-genital Contact TotalLen": 7000, "Side by side Contact TotalLen": 10000,
               "Side by side Contact, opposite way TotalLen": 10000, "Train2 TotalLen": 200,
               "FollowZone Isolated TotalLen": 1500, "Social approach TotalLen": 10000,
               "Approach contact TotalLen": 10000, "Get away TotalLen": 12000, "Break contact TotalLen": 3000,

               "Move isolated MeanDur": 15, "Move in contact MeanDur": 10, "WallJump MeanDur": 4,
               "Stop isolated MeanDur": 180, "Rear isolated MeanDur": 12,
               "Rear in contact MeanDur": 12, "Contact MeanDur": 150, "Group2 MeanDur": 200,
               "Group3 MeanDur": 200, "Oral-oral Contact MeanDur": 10,
               "Oral-genital Contact MeanDur": 10, "Side by side Contact MeanDur": 15,
               "Side by side Contact, opposite way MeanDur": 15, "Train2 MeanDur": 20,
               "FollowZone Isolated MeanDur": 6, "Social approach MeanDur": 6,
               "Approach contact MeanDur": 6, "Get away MeanDur": 6, "Break contact MeanDur": 30,

               "Move isolated Nb": 80000, "Move in contact Nb": 40000, "WallJump Nb": 6000,
               "Stop isolated Nb": 80000, "Rear isolated Nb": 50000,
               "Rear in contact Nb": 40000, "Contact Nb": 25000, "Group2 Nb": 100000,
               "Group3 Nb": 100000, "Oral-oral Contact Nb": 40000,
               "Oral-genital Contact Nb": 40000, "Side by side Contact Nb": 40000,
               "Side by side Contact, opposite way Nb": 40000, "Train2 Nb": 600,
               "FollowZone Isolated Nb": 25000, "Social approach Nb": 140000,
               "Approach contact Nb": 100000, "Get away Nb": 100000, "Break contact Nb": 8000}

    xPos = [0, 1, 2]
    yLabelText = {' TotalLen': ' duration (s)', ' MeanDur': ' mean duration (frames)', ' Nb': ' occurrences'}
    letterList = list(string.ascii_uppercase)

    while True:

        question = "Do you want to:"
        question += "\n\t [1] plot the behavioral profiles for C57BL/6J mice (males and females, at the different ages)?"
        question += "\n\t [2] plot the behavioral profiles for Shank3 mice (females)?"
        question += "\n\t [3] plot the figure 2 with the behavioral profiles for C57BL/6J mice and for Shank3 mice?"
        question += "\n"
        answer = inpuinputFilestion)

        if answer == "1":
            #plot the complete figure for the profiles
            #open the json file with behavioral profile data for pairs
            jsonFile = getJsonFileToProcess()

            with open(jsonFile) as json_data:
                dataByNight = json.load(json_data)
            print("json file re-imported.")

            #Merge profile data over the three nights
            data = mergeProfileOverNights(profileData=dataByNight, categoryList=[' TotalLen', ' Nb', ' MeanDur'], behaviouralEventOneMouse=behaviouralEventOneMouse)

            #Create a data dictionary with all the behaviors of interest
            dataDic = {}

            eventSingleList = ['totalDistance', "Move isolated MeanDur", "Move in contact MeanDur", "WallJump MeanDur",
            "Stop isolated MeanDur", "Rear isolated MeanDur",
            "Rear in contact MeanDur", "Oral-genital Contact MeanDur", "Train2 MeanDur", "FollowZone Isolated MeanDur", "Social approach MeanDur",
            "Approach contact MeanDur", "Get away MeanDur", "Break contact MeanDur",
            "Move isolated Nb", "Move in contact Nb", "WallJump Nb",
            "Stop isolated Nb", "Rear isolated Nb",
            "Rear in contact Nb", "Oral-genital Contact Nb", "Train2 Nb", "FollowZone Isolated Nb", "Social approach Nb",
            "Approach contact Nb", "Get away Nb", "Break contact Nb",
            "Move isolated TotalLen", "Move in contact TotalLen", "WallJump TotalLen",
            "Stop isolated TotalLen", "Rear isolated TotalLen",
            "Rear in contact TotalLen", "Oral-genital Contact TotalLen", "Train2 TotalLen", "FollowZone Isolated TotalLen", "Social approach TotalLen",
            "Approach contact TotalLen", "Get away TotalLen", "Break contact TotalLen"]

            for event in eventSingleList:
                print('##Event: ', event)
                dataDic[event] = getProfileValues(profileData=data, night='all nights', event=event)
                #print(dataDic[event])

            eventPairList = ["Contact MeanDur", "Oral-oral Contact MeanDur",  "Side by side Contact MeanDur", "Side by side Contact, opposite way MeanDur",
            "Contact Nb", "Oral-oral Contact Nb", "Side by side Contact Nb", "Side by side Contact, opposite way Nb",
            "Contact TotalLen", "Oral-oral Contact TotalLen", "Side by side Contact TotalLen", "Side by side Contact, opposite way TotalLen"]

            for event in eventPairList:
                print('##Event: ', event)
                dataDic[event] = getProfileValuesPairs(profileData=data, night='all nights', event=event)
                #print(dataDic[event])


            for valueCat in [' TotalLen', ' Nb', ' MeanDur']:

                # create the figure
                nbCol = 4
                nbRow = 4
                row = 0
                col = 0
                fig, axes = plt.subplots(nrows=nbRow, ncols=nbCol, figsize=(4*nbCol, 3*nbRow), sharey=False)
                k = 0

                ##############################
                for behavEvent in behaviouralEventOneMouseShort:
                    ax = axes[row, col]
                    # create a dataframe from the dictionary to plot data and test them
                    event = behavEvent + valueCat
                    df = pd.DataFrame.from_dict(dataDic[event])
                    # plot data
                    singlePlotPerEventProfileBothSexesPerAge(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                             ax=ax, letter=letterList[k],
                                                             yMin=yMinDic[behavEvent + valueCat],
                                                             yMax=yMaxDic[behavEvent + valueCat], xPos=xPos)
                    ax.set_ylabel(yLabelText[valueCat], fontsize=12)
                    ax.set_title('{}'.format(getFigureBehaviouralEventsLabels(behavEvent)), fontsize=13, pad=0.2, fontweight='bold', position=(0.5,1.6))
                    ax.get_xaxis().get_label().set_visible(False)

                    # test differences between sexes for each age class:
                    testProfileDataBetweenSexesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic, yMaxDic=yMaxDic)
                    # test differences between age classes within sex:
                    testProfileDataBetweenAgesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic, yMaxDic=yMaxDic)
                    k += 1
                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1
                valCatGeneral = valueCat

                # for distance, plot the data
                if valCatGeneral == ' TotalLen':
                    ax = axes[row, col]
                    behavEvent = 'totalDistance'
                    valueCat = ''
                    # create a dataframe from the dictionary to plot data and test them
                    event = behavEvent + valueCat
                    df = pd.DataFrame.from_dict(dataDic[event])
                    # plot data
                    singlePlotPerEventProfileBothSexesPerAge(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                             ax=ax, letter=letterList[k],
                                                             yMin=yMinDic[behavEvent + valueCat],
                                                             yMax=yMaxDic[behavEvent + valueCat], xPos=xPos)
                    ax.set_ylabel('distance travelled (m)', fontsize=12)
                    ax.set_title('distance travelled', fontsize=13, pad=0.2,
                                 fontweight='bold', position=(0.5,1.6))
                    ax.get_xaxis().get_label().set_visible(False)

                    # test differences between sexes for each age class:
                    testProfileDataBetweenSexesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                                  yMaxDic=yMaxDic)
                    # test differences between age classes within sex:
                    testProfileDataBetweenAgesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                                 yMaxDic=yMaxDic)
                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1
                    ##############################
                    ax = axes[row, col]
                    ax.axis('off')

                    malePatch = {}
                    femalePatch = {}
                    for age in ageList:
                        malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age),
                                                        label='male {}'.format(age))
                        femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age),
                                                          label='female {}'.format(age))

                    handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'], malePatch['7mo'],
                               femalePatch['7mo']]
                    ax.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1

                    ################################
                    ax = axes[row, col]
                    ax.axis('off')

                else:
                    ax = axes[row, col]
                    ax.axis('off')

                    malePatch = {}
                    femalePatch = {}
                    for age in ageList:
                        malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age),
                                                        label='male {}'.format(age))
                        femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age),
                                                          label='female {}'.format(age))

                    handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'],
                               malePatch['7mo'],
                               femalePatch['7mo']]
                    ax.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1

                    ################################
                    ax = axes[row, col]
                    ax.axis('off')

                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1

                    ################################
                    ax = axes[row, col]
                    ax.axis('off')

                fig.show()
                fig.tight_layout()
                fig.savefig('profile_complete_pairs_age_sex_{}.pdf'.format(valCatGeneral), dpi=300)
                fig.savefig('profile_complete_pairs_age_sex_{}.jpg'.format(valCatGeneral), dpi=300)

            print('Job done.')
            break


        if answer == '2':
            # plot the complete figure for the profiles
            # open the json file with behavioral profile data for pairs
            jsonFile = getJsonFileToProcess()

            with open(jsonFile) as json_data:
                dataByNight = json.load(json_data)
            print("json file re-imported.")

            #Merge profile data over the three nights
            data = mergeProfileOverNights(profileData=dataByNight, categoryList=[' TotalLen', ' Nb', ' MeanDur'],
                                          behaviouralEventOneMouse=behaviouralEventOneMouse)

            #data = dataByNight
            # Create a data dictionary with all the behaviors of interest
            dataDic = {}

            eventSingleList = ['totalDistance', "Move isolated MeanDur", "Move in contact MeanDur", "WallJump MeanDur",
                               "Stop isolated MeanDur", "Rear isolated MeanDur",
                               "Rear in contact MeanDur", "Oral-genital Contact MeanDur", "Train2 MeanDur",
                               "FollowZone Isolated MeanDur", "Social approach MeanDur",
                               "Approach contact MeanDur", "Get away MeanDur", "Break contact MeanDur",
                               "Move isolated Nb", "Move in contact Nb", "WallJump Nb",
                               "Stop isolated Nb", "Rear isolated Nb",
                               "Rear in contact Nb", "Oral-genital Contact Nb", "Train2 Nb", "FollowZone Isolated Nb",
                               "Social approach Nb",
                               "Approach contact Nb", "Get away Nb", "Break contact Nb",
                               "Move isolated TotalLen", "Move in contact TotalLen", "WallJump TotalLen",
                               "Stop isolated TotalLen", "Rear isolated TotalLen",
                               "Rear in contact TotalLen", "Oral-genital Contact TotalLen", "Train2 TotalLen",
                               "FollowZone Isolated TotalLen", "Social approach TotalLen",
                               "Approach contact TotalLen", "Get away TotalLen", "Break contact TotalLen"]

            for event in eventSingleList:
                print('##Event: ', event)
                dataDic[event] = getProfileValues(profileData=data, night='all nights', event=event)
                # print(dataDic[event])

            eventPairList = ["Contact MeanDur", "Oral-oral Contact MeanDur", "Side by side Contact MeanDur",
                             "Side by side Contact, opposite way MeanDur",
                             "Contact Nb", "Oral-oral Contact Nb", "Side by side Contact Nb",
                             "Side by side Contact, opposite way Nb",
                             "Contact TotalLen", "Oral-oral Contact TotalLen", "Side by side Contact TotalLen",
                             "Side by side Contact, opposite way TotalLen"]

            for event in eventPairList:
                print('##Event: ', event)
                dataDic[event] = getProfileValuesPairs(profileData=data, night='all nights', event=event)
                # print(dataDic[event])

            for valueCat in [' TotalLen', ' Nb', ' MeanDur']:

                # create the figure
                nbCol = 5
                nbRow = 3
                row = 0
                col = 0
                fig, axes = plt.subplots(nrows=nbRow, ncols=nbCol, figsize=(3 * nbCol, 3 * nbRow), sharey=False)
                k = 0

                ##############################
                for behavEvent in behaviouralEventOneMouseShort:
                    ax = axes[row, col]
                    # create a dataframe from the dictionary to plot data and test them
                    event = behavEvent + valueCat
                    df = pd.DataFrame.from_dict(dataDic[event])
                    # plot data
                    singlePlotPerEventProfilePerGenotype(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                             ax=ax, letter=letterList[k],
                                                             yMin=yMinDic[behavEvent + valueCat],
                                                             yMax=yMaxDic[behavEvent + valueCat], xPos=[0])
                    ax.set_ylabel(yLabelText[valueCat], fontsize=12)
                    ax.set_title('{}'.format(getFigureBehaviouralEventsLabels(behavEvent)), fontsize=13,
                                 fontweight='bold', position=(0.5,1.6))
                    ax.get_xaxis().get_label().set_visible(False)

                    # test differences between sexes for each age class:
                    testProfileDataBetweenGenotypes(df=df, ax=ax, event=event, yMinDic=yMinDic,
                                                  yMaxDic=yMaxDic)

                    k += 1
                    if col < 4:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1
                valCatGeneral = valueCat

                # for distance, plot the data
                if valCatGeneral == ' TotalLen':
                    ax = axes[row, col]
                    behavEvent = 'totalDistance'
                    valueCat = ''
                    # create a dataframe from the dictionary to plot data and test them
                    event = behavEvent + valueCat
                    df = pd.DataFrame.from_dict(dataDic[event])
                    # plot data
                    singlePlotPerEventProfilePerGenotype(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                             ax=ax, letter=letterList[k],
                                                             yMin=yMinDic[behavEvent + valueCat],
                                                             yMax=yMaxDic[behavEvent + valueCat], xPos=[0])
                    ax.set_ylabel('distance travelled (m)', fontsize=12)
                    ax.set_title('distance travelled', fontsize=13,
                                 fontweight='bold', position=(0.5,1.6))
                    ax.get_xaxis().get_label().set_visible(False)

                    # test differences between sexes for each age class:
                    testProfileDataBetweenGenotypes(df=df, ax=ax, event=event, yMinDic=yMinDic,
                                                  yMaxDic=yMaxDic)

                    if col < 4:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1
                    ##############################
                    ax = axes[row, col]
                    ax.axis('off')

                    wtPatch = mpatches.Patch(edgecolor='black', facecolor=getColorWT(),
                                                        label='C57BL/6J female')
                    koPatch = mpatches.Patch(edgecolor='black', facecolor=getColorKO(),
                                                          label='Shank3-/- female')

                    handles = [wtPatch, koPatch]
                    ax.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

                    if col < 3:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1


                else:
                    ax = axes[row, col]
                    ax.axis('off')

                    wtPatch = mpatches.Patch(edgecolor='black', facecolor=getColorWT(),
                                             label='C57BL/6J female')
                    koPatch = mpatches.Patch(edgecolor='black', facecolor=getColorKO(),
                                             label='Shank3-/- female')

                    handles = [wtPatch, koPatch]
                    ax.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

                    if col < 4:
                        col += 1
                        row = row
                    else:
                        col = 0
                        row += 1

                    ################################
                    ax = axes[row, col]
                    ax.axis('off')


                fig.show()
                fig.tight_layout()
                fig.savefig('profile_complete_pairs_B6_shank3_{}.pdf'.format(valCatGeneral), dpi=300)
                fig.savefig('profile_complete_pairs_B6_shank3_{}.jpg'.format(valCatGeneral), dpi=300)

            print('Job done.')
            break

        if answer == "3":
            # plot the complete figure for the profiles
            # open the json file with behavioral profile data for pairs
            jsonFile = getJsonFileToProcess()

            with open(jsonFile) as json_data:
                dataByNight = json.load(json_data)
            print("json file re-imported.")

            # Merge profile data over the three nights
            data = mergeProfileOverNights(profileData=dataByNight, categoryList=[' TotalLen', ' Nb', ' MeanDur'],
                                          behaviouralEventOneMouse=behaviouralEventOneMouse)

            # Create a data dictionary with all the behaviors of interest
            dataDic = {}

            eventSingleList = ['totalDistance', "Move isolated MeanDur", "Move in contact MeanDur", "WallJump MeanDur",
                               "Stop isolated MeanDur", "Rear isolated MeanDur",
                               "Rear in contact MeanDur", "Oral-genital Contact MeanDur", "Train2 MeanDur",
                               "FollowZone Isolated MeanDur", "Social approach MeanDur",
                               "Approach contact MeanDur", "Get away MeanDur", "Break contact MeanDur",
                               "Move isolated Nb", "Move in contact Nb", "WallJump Nb",
                               "Stop isolated Nb", "Rear isolated Nb",
                               "Rear in contact Nb", "Oral-genital Contact Nb", "Train2 Nb", "FollowZone Isolated Nb",
                               "Social approach Nb",
                               "Approach contact Nb", "Get away Nb", "Break contact Nb",
                               "Move isolated TotalLen", "Move in contact TotalLen", "WallJump TotalLen",
                               "Stop isolated TotalLen", "Rear isolated TotalLen",
                               "Rear in contact TotalLen", "Oral-genital Contact TotalLen", "Train2 TotalLen",
                               "FollowZone Isolated TotalLen", "Social approach TotalLen",
                               "Approach contact TotalLen", "Get away TotalLen", "Break contact TotalLen"]

            for event in eventSingleList:
                print('##Event: ', event)
                dataDic[event] = getProfileValues(profileData=data, night='all nights', event=event)
                # print(dataDic[event])

            eventPairList = ["Contact MeanDur", "Oral-oral Contact MeanDur", "Side by side Contact MeanDur",
                             "Side by side Contact, opposite way MeanDur",
                             "Contact Nb", "Oral-oral Contact Nb", "Side by side Contact Nb",
                             "Side by side Contact, opposite way Nb",
                             "Contact TotalLen", "Oral-oral Contact TotalLen", "Side by side Contact TotalLen",
                             "Side by side Contact, opposite way TotalLen"]

            for event in eventPairList:
                print('##Event: ', event)
                dataDic[event] = getProfileValuesPairs(profileData=data, night='all nights', event=event)
                # print(dataDic[event])

            # create the figure
            nbCol = 4
            nbRow = 3
            row = 0
            col = 0
            fig, axes = plt.subplots(nrows=nbRow, ncols=nbCol, figsize=(3.5 * nbCol, 3 * nbRow), sharey=False)
            k = 0

            behaviourListWt = ['Oral-genital Contact', 'Contact', 'Oral-oral Contact',
                            'Oral-genital Contact', 'Approach contact', 'Get away']
            valueCatListWt = [' TotalLen', ' MeanDur', ' MeanDur',
                            ' MeanDur', ' MeanDur', ' MeanDur']

            behaviourListKo = ['Contact', 'Oral-oral Contact', 'Side by side Contact']
            valueCatListKo = [' TotalLen', ' MeanDur', ' MeanDur']

            ax = axes[row, col]
            behavEvent = 'totalDistance'
            valueCat = ''
            # create a dataframe from the dictionary to plot data and test them
            event = behavEvent + valueCat
            df = pd.DataFrame.from_dict(dataDic[event])
            # plot data
            singlePlotPerEventProfileBothSexesPerAge(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                     ax=ax, letter=letterList[k],
                                                     yMin=yMinDic[behavEvent + valueCat],
                                                     yMax=yMaxDic[behavEvent + valueCat], xPos=xPos)
            ax.set_ylabel('distance travelled (m)', fontsize=12)
            ax.set_title('distance travelled', fontsize=12, pad=0.2,
                         fontweight='bold')
            ax.get_xaxis().get_label().set_visible(False)
            ax.text(x=-0.3, y=yMinDic[behavEvent + valueCat], s='M')
            ax.text(x=0.2, y=yMinDic[behavEvent + valueCat], s='F')
            ax.text(x=0.7, y=yMinDic[behavEvent + valueCat], s='M')
            ax.text(x=1.2, y=yMinDic[behavEvent + valueCat], s='F')
            ax.text(x=1.7, y=yMinDic[behavEvent + valueCat], s='M')
            ax.text(x=2.2, y=yMinDic[behavEvent + valueCat], s='F')

            # test differences between sexes for each age class:
            testProfileDataBetweenSexesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                          yMaxDic=yMaxDic)
            # test differences between age classes within sex:
            testProfileDataBetweenAgesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                         yMaxDic=yMaxDic)
            k += 1
            if col < 3:
                col += 1
                row = row
            else:
                col = 0
                row += 1

            ##############################
            b = 0
            for behavEvent in behaviourListWt:
                ax = axes[row, col]
                # create a dataframe from the dictionary to plot data and test them
                valueCat = valueCatListWt[b]
                event = behavEvent + valueCat
                df = pd.DataFrame.from_dict(dataDic[event])
                # plot data
                singlePlotPerEventProfileBothSexesPerAge(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                         ax=ax, letter=letterList[k],
                                                         yMin=yMinDic[behavEvent + valueCat],
                                                         yMax=yMaxDic[behavEvent + valueCat], xPos=xPos)
                ax.set_ylabel(yLabelText[valueCat], fontsize=12)
                ax.set_title('{}'.format(getFigureBehaviouralEventsLabels(behavEvent)), fontsize=12, pad=0.2,
                             fontweight='bold')
                ax.get_xaxis().get_label().set_visible(False)

                # test differences between sexes for each age class:
                testProfileDataBetweenSexesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                              yMaxDic=yMaxDic)
                # test differences between age classes within sex:
                testProfileDataBetweenAgesWT(ageList=ageList, df=df, ax=ax, event=event, yMinDic=yMinDic,
                                             yMaxDic=yMaxDic)
                k += 1
                b += 1
                if col < 3:
                    col += 1
                    row = row
                else:
                    col = 0
                    row += 1

            ##################################
            ax = axes[row, col]
            ax.axis('off')
            malePatch = {}
            femalePatch = {}
            for age in ageList:
                malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age), linewidth=1.2,
                                                label='B6 male {}'.format(age))
                femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age), linewidth=1.2,
                                                  label='B6 female {}'.format(age))
            koPatch = mpatches.Patch(edgecolor='black', facecolor=getColorKO(), linewidth=1.2,
                                     label='Shank3-/- female 3mo')

            handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'],
                       malePatch['7mo'], femalePatch['7mo'], koPatch]
            ax.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

            if col < 3:
                col += 1
                row = row
            else:
                col = 0
                row += 1


            ##############################
            b = 0
            for behavEvent in behaviourListKo:
                ax = axes[row, col]
                # create a dataframe from the dictionary to plot data and test them
                valueCat = valueCatListKo[b]
                event = behavEvent + valueCat
                df = pd.DataFrame.from_dict(dataDic[event])
                # plot data
                singlePlotPerEventProfilePerGenotype(dataframe=df, valueCat=valueCat, behavEvent=behavEvent,
                                                         ax=ax, letter=letterList[k],
                                                         yMin=yMinDic[behavEvent + valueCat],
                                                         yMax=yMaxDic[behavEvent + valueCat], xPos=[0])
                ax.set_ylabel(yLabelText[valueCat], fontsize=12)
                ax.set_title('{}'.format(getFigureBehaviouralEventsLabels(behavEvent)), fontsize=12, pad=0.2,
                             fontweight='bold')
                ax.get_xaxis().get_label().set_visible(False)
                if b == 0:
                    ax.text(x=-0.16, y=yMinDic[behavEvent + valueCat]+0.04 * (yMaxDic[behavEvent + valueCat] - yMinDic[behavEvent + valueCat]), s='B6 F')
                    ax.text(x=0.07, y=yMinDic[behavEvent + valueCat]+0.04 * (yMaxDic[behavEvent + valueCat] - yMinDic[behavEvent + valueCat]), s='Shank3 F')

                # test differences between genotype:
                testProfileDataBetweenGenotypes(df=df, ax=ax, event=event, yMinDic=yMinDic,
                                              yMaxDic=yMaxDic)

                k += 1
                b += 1
                if col < 3:
                    col += 1
                    row = row
                else:
                    col = 0
                    row += 1

            ##################
            ax = axes[row, col]
            ax.axis('off')

            fig.show()
            fig.tight_layout()
            fig.savefig('Figure2_behavioral_profile_complete_pairs_B6_shank3.pdf', dpi=300)
            fig.savefig('Figure2_behavioral_profile_complete_pairs_B6_shank3.jpg', dpi=300)

            print('Job done.')
            break

