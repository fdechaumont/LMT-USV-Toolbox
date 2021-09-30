'''
Created on 20 janv. 2020

@author: Elodie
'''

from tkinter.filedialog import askopenfilename
import sqlite3
import os
import numpy as np
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
from scipy.stats import mannwhitneyu, kruskal, ttest_ind
from numpy import sum
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches

import numpy as np

from scipy import signal
from scipy.io import wavfile
import os
import wave
import pylab
from lmtanalysis.FileUtil import getFilesToProcess
import pandas as pd
from pandas.core.frame import DataFrame
from scipy import stats

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import json
from LMT.USV.testVocTraitUsagePerEventType.test1 import cleanVoc
from LMT.USV.lib.vocUtil import *
from LMT.USV.experimentList.experimentList import getAllExperimentList
from LMT.USV.figure.figParameter import *
from USV.figure.figUtil import addJitter, getStarsFromPvalues
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import string
from LMT.USV.testVocTraitUsagePerEventType.AcousticAnalysesAllUsv import *


if __name__ == '__main__':
    '''
    This codes plots the figure panels about the acoustic variables of USV. Data are extracted from a json file.
    '''
    pd.set_option('display.max_rows', 20)
    pd.set_option('display.max_columns', 20)
    #pd.set_option('display.width', 1000)

    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})
    acousticVariables = getFigureVocTraits()
    colNames = ['strain', 'sex', 'age', 'pair'] + acousticVariables
    colSex = {'female': 'black', 'male': 'grey'}
    print(colNames)

    letterList = list(string.ascii_uppercase)
    k= 0

    acousticVariablesSelectedWt = ['durationMs', 'nbModulation', 'griffIndex']
    acousticVariablesSelectedKo = ['durationMs', 'nbModulation', 'nbJump', 'startFrequencyHz', 'meanPower']

    with open("dataAcousticAnalysisAllUsvs.json") as json_data:
        dataDf = json.load(json_data)
    print("json file for acoustic variables re-imported.")
    df = pd.DataFrame(dataDf, columns=colNames)
    print(df.head(n=10))
    print('Dataframe created.')

    # Create general figure
    gs = gridspec.GridSpec(3, 15)
    fig = plt.figure(figsize=(18, 9))

    # spectro of male USVs 3mo
    letter = letterList[k]
    ax = fig.add_subplot(gs[0, 0:5])
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.patch.set_facecolor('white')

    imgPos = (0.4, 0.5)
    image = 'spectro_m_3mo.png'
    behavSchema = mpimg.imread(image)
    imgBox = OffsetImage(behavSchema, zoom=0.27)
    imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
    ax.add_artist(imageBox)

    ax.text(-0.15, 1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    ax.set_title(label='B6 males (3mo)', loc='center', fontsize=12, horizontalalignment='center', color='black', weight='bold')

    k += 1

    # spectro of female USVs 3mo
    letter = letterList[k]
    ax = fig.add_subplot(gs[0, 5:10])
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_facecolor('white')

    imgPos = (0.35, 0.5)
    image = 'spectro_f_3mo.png'
    behavSchema = mpimg.imread(image)
    imgBox = OffsetImage(behavSchema, zoom=0.27)
    imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
    ax.add_artist(imageBox)

    ax.text(-0.18, 1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    ax.set_title(label='B6 females (3mo)', loc='center', fontsize=12, horizontalalignment='center', color='black', weight='bold')
    k += 1

    # spectro of female USVs 3mo
    letter = letterList[k]
    ax = fig.add_subplot(gs[0, 10:15])
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_facecolor('white')

    imgPos = (0.45, 0.5)
    image = 'spectro_s3_3mo.png'
    behavSchema = mpimg.imread(image)
    imgBox = OffsetImage(behavSchema, zoom=0.27)
    imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
    ax.add_artist(imageBox)

    ax.text(-0.2, 1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')
    ax.set_title(label='Shank3-/- females (3mo)', loc='center', fontsize=12, horizontalalignment='center', color='black', weight='bold')
    k += 1

    ################################################
    #boxplots for B6 data, with both sexes and all age classes for raw values
    # create a panda dataframe with acoustic measures for all voc
    xPos = {'5we': 0, '3mo': 1, '7mo': 2}
    # select the data only for C57BL/6J mice
    dfWT = df.loc[df['strain'] == 'C57BL/6J', :]

    print('################################')
    bp = plotAcousticVarWTBoxplot(ax=fig.add_subplot(gs[1, 0:5]), dataframe=dfWT, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=acousticVariablesSelectedWt[0], letter=letterList[k])
    k += 1
    malePatch = {}
    femalePatch = {}
    for age in ageList:
        malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age), label='male {}'.format(age))
        femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age), label='female {}'.format(age))

    handles=[malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'], malePatch['7mo'], femalePatch['7mo']]
    bp.legend(handles=handles, loc=(0.05, 0.3)).set_visible(True)

    print('################################')
    plotAcousticVarWTBoxplot(ax=fig.add_subplot(gs[1, 5:10]), dataframe=dfWT, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedWt[1], letter=letterList[k])
    k += 1

    print('################################')
    plotAcousticVarWTBoxplot(ax=fig.add_subplot(gs[1, 10:15]), dataframe=dfWT, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedWt[2], letter=letterList[k])
    k += 1


    '''################################################
    # age effects for B6 data, with both sexes for raw values
    colSex = {'male': 'grey', 'female': 'black'}
    print('################################')
    plotAcousticVarAgeEffect(ax=fig.add_subplot(gs[2, 0:5]), dfWT=dfWT, var=acousticVariablesSelectedWt[0], yMinVarAge=yMinVarAge, yMaxVarAge=yMaxVarAge, ageList=ageList,
                             colSex=colSex, letter=letterList[k])
    malePatch = mpatches.Patch(edgecolor=colSex['male'], facecolor='white', label='male')
    femalePatch = mpatches.Patch(edgecolor=colSex['female'], facecolor='white', label='female')
    handles = [malePatch, femalePatch]
    plt.legend(handles=handles, loc=(0.05, 0.6)).set_visible(True)
    k += 1

    print('################################')
    plotAcousticVarAgeEffect(ax=fig.add_subplot(gs[2, 5:10]), dfWT=dfWT, var=acousticVariablesSelectedWt[1],
                             yMinVarAge=yMinVarAge, yMaxVarAge=yMaxVarAge, ageList=ageList,
                             colSex=colSex, letter=letterList[k])
    k += 1

    print('################################')
    plotAcousticVarAgeEffect(ax=fig.add_subplot(gs[2, 10:15]), dfWT=dfWT, var=acousticVariablesSelectedWt[2],
                             yMinVarAge=yMinVarAge, yMaxVarAge=yMaxVarAge, ageList=ageList,
                             colSex=colSex, letter=letterList[k])
    k += 1
'''
    ################################################
    # boxplots for genotype effects at 3 mo for raw values in females
    # select the data only for C57BL/6J and Shank3 3mo female mice
    dfKO = df.loc[(df['sex'] == 'female') & (df['age'] == '3mo'), :]

    print('################################')
    bp = plotAcousticVarKOBoxplot(ax=fig.add_subplot(gs[2, 0:3]), dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=acousticVariablesSelectedKo[0],
                             letter=letterList[k])
    wtPatch = mpatches.Patch(edgecolor=colSex['female'], facecolor=getColorWT(), label='B6 female')
    koPatch = mpatches.Patch(edgecolor=colSex['female'], facecolor=getColorKO(), label='Shank3 female')
    handles = [wtPatch, koPatch]
    bp.legend(handles=handles, loc=(0.05, 0.7)).set_visible(True)
    k += 1

    print('################################')
    plotAcousticVarKOBoxplot(ax=fig.add_subplot(gs[2, 3:6]), dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedKo[1], letter=letterList[k])
    k += 1

    print('################################')
    plotAcousticVarKOBoxplot(ax=fig.add_subplot(gs[2, 6:9]), dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedKo[2], letter=letterList[k])
    k += 1

    print('################################')
    plotAcousticVarKOBoxplot(ax=fig.add_subplot(gs[2, 9:12]), dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedKo[3], letter=letterList[k])
    k += 1

    print('################################')
    plotAcousticVarKOBoxplot(ax=fig.add_subplot(gs[2, 12:15]), dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar,
                             variable=acousticVariablesSelectedKo[4], letter=letterList[k])
    k += 1



    plt.show()
    gs.tight_layout(fig)
    fig.savefig('figure4_acoustic_var.pdf', dpi=300)
    fig.savefig('figure4_acoustic_var.jpg', dpi=300)

    print('Job done.')



