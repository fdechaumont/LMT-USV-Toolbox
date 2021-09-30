'''
Created by E. Ey on 07/12/2020
'''

from lmtanalysis.Event import *
from lmtanalysis.Measure import *
import numpy as np; np.random.seed(0)
from tkinter.filedialog import askopenfilename
from lmtanalysis.Util import getMinTMaxTAndFileNameInput, getMinTMaxTInput
import sqlite3
from lmtanalysis.FileUtil import getFilesToProcess
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
from LMT.USV.figure.figUtil import addJitter
from LMT.USV.figure.figParameter import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, kruskal, ttest_ind
import matplotlib.image as mpimg
import statsmodels.formula.api as smf

from LMT.USV.usvDescription.Compute_Number_USVs import *
from LMT.USV.usvEventsCorrelations.Compute_Correlation_USV_Burst_With_Events import *
from LMT.USV.figure.burstTraitUsagePerEventContext import *
from LMT.USV.usvDescription.Compute_Number_USVs import computeNumberUsv,\
    createDataframeFromJson, getDataFrameWT, getDataFrameKO,\
    plotNumberUsvDayNight, plotNumberUsvWTWithAge, plotNumberUsvKO,\
    createDataframeFromJsonNumberUsvPerBurst, plotNumberUsvPerBurstWTBoxplot,\
    plotNumberUsvPerBurstKOBoxplot
from LMT.USV.figure.figParameter import getFigureBehaviouralEvents, getColorAge,\
    getColorWT, getColorKO
from LMT.USV.usvEventsCorrelations.Compute_Correlation_USV_Burst_With_Events import generateDataframeCorrelationFromDic,\
    plotCorrelationUsvEvents, plotCorrelationUsvEventsWithUSVWithKo
from LMT.USV.figure.burstTraitUsagePerEventContext import plotBurstTraitsBoxplotsPerAgePerContextsPerGeno,\
    plotBurstTraitUsagePerEventContextPerSet2
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches

if __name__ == '__main__':

    print("Code launched.")

    # set font
    from matplotlib import rc, gridspec

    rc('font', **{'family': 'serif', 'serif': ['Arial']})

    while True:

        question = "Do you want to:"
        question += "\n\t [1] compute the total number of USVs, the number of USV bursts and the number of USV per burst over the three days?"
        question += "\n\t [2] plot figure 3 with data for C57BL/6J mice and for Shank3 mice?"
        question += "\n\t [3] plot figure 7 with bursts characteristics with data for C57BL/6J mice and for Shank3 mice?"
        question += "\n\t [4] plot suppl figure with heatmaps for bursts characteristics with data for C57BL/6J mice and for Shank3 mice?"
        question += "\n"
        answer = inputFile(question)

        if answer == "1":
            computeNumberUsv(tmin=0, tmax=7776000)
            print("ok")

            break

        if answer == "2":
            strainList = ['C57BL/6J', 'Shank3']
            ageList = ['5we', '3mo', '7mo']
            sexList = ['male', 'female']
            eventListToTest = getFigureBehaviouralEvents(longList=True)
            df = createDataframeFromJson(jsonFile='dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList)
            # print(df)
            ###################################################
            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')

            # data KO
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')

            #################################################
            for var in ['nbUsv', 'nbUsvDay', 'nbUsvNight', 'nbBurst', 'nbUsvBurst']:
                print('tests for age effect: ')
                for sexClass in sexList:
                    K, p = stats.kruskal(
                        dataframeWT[var][(dataframeWT['age'] == '5we') & (dataframeWT['sex'] == sexClass)],
                        dataframeWT[var][(dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == sexClass)],
                        dataframeWT[var][(dataframeWT['age'] == '7mo') & (dataframeWT['sex'] == sexClass)])
                    print(var, ' ', sexClass, ' K = ', K, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age'] == '5we') & (dataframeWT['sex'] == sexClass)],
                                    dataframeWT[var][(dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == sexClass)],
                                    alternative='two-sided')
                    print(var, ' ', sexClass, ' 5we vs 3mo', ' U = ', U, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age'] == '3mo') & (dataframeWT['sex'] == sexClass)],
                                    dataframeWT[var][(dataframeWT['age'] == '7mo') & (dataframeWT['sex'] == sexClass)],
                                    alternative='two-sided')
                    print(var, ' ', sexClass, ' 3mo vs 7mo', ' U = ', U, 'p = ', p)
                    U, p = wilcoxon(dataframeWT[var][(dataframeWT['age'] == '5we') & (dataframeWT['sex'] == sexClass)],
                                    dataframeWT[var][(dataframeWT['age'] == '7mo') & (dataframeWT['sex'] == sexClass)],
                                    alternative='two-sided')
                    print(var, ' ', sexClass, ' 5we vs 7mo', ' U = ', U, 'p = ', p)

                print('Tests between sexes:')
                for ageClass in ageList:
                    U, p = mannwhitneyu(
                        dataframeWT[var][(dataframeWT['age'] == ageClass) & (dataframeWT['sex'] == 'male')],
                        dataframeWT[var][(dataframeWT['age'] == ageClass) & (dataframeWT['sex'] == 'female')],
                        alternative='two-sided')
                    print(var, ' ', ageClass, ' U = ', U, 'p = ', p)

            #########################################################
            # plots
            yMinUsv = {'nbUsv': 0, 'nbBurst': 0, 'nbUsvBurst': 0}
            yMaxUsv = {'nbUsv': 24000, 'nbBurst': 5000, 'nbUsvBurst': 50}
            jitterValue = 0.2
            x1 = {'nbUsvDay': {'5we': 1, '3mo': 4, '7mo': 7}, 'nbUsvNight': {'5we': 2, '3mo': 5, '7mo': 8},
                  'nbBurst': {'5we': 1.5, '3mo': 4.5, '7mo': 7.5}, 'nbUsvBurst': {'5we': 1.5, '3mo': 4.5, '7mo': 7.5}}

            # Create general figure
            gs = gridspec.GridSpec(2, 29)
            # fig = plt.subplots(figsize=(24, 15), sharex=False, sharey=False)
            fig = plt.figure(figsize=(16, 3))

            ################################################
            # plot the number of usv in days and nights
            df = createDataframeFromJson(jsonFile='dataUsvDescription750.json', strainList=strainList, ageList=ageList, sexList=sexList)
            # dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            plotNumberUsvDayNight(ax=fig.add_subplot(gs[:, 0:3]), dataframeWT=df, yMinUsv=yMinUsv,
                                  yMaxUsv=yMaxUsv, jitterValue=jitterValue, letter='A', sexList=sexList, ageList=ageList)
            print('day night done')

            ################################################
            # plot the number of usv in males versus females with all ages
            plotNumberUsvWTWithAge(ax=fig.add_subplot(gs[:, 5:9]), dataframeWT=dataframeWT, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, jitterValue=jitterValue, letter='B', ageList=ageList)
            #plotNumberUsvWT(ax=fig.add_subplot(gs[0, 0:1]), dataframeWT=dataframeWT, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, jitterValue=jitterValue, letter='A')

            ################################################
            # plot the number of usv in C57BL/6J versus Shank3 KO females
            plotNumberUsvKO(ax=fig.add_subplot(gs[:, 11:13]), dataframeKO=dataframeKO, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, jitterValue=jitterValue, letter='C')


            ################################################
            # plot the correlation of usv with events
            # open the json file to work on pre-computed data
            jsonFile = "dataCorrelationRefUsvBurst.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file re-imported.")
            selectedEventList = ['Stop isolated', 'Contact', 'Approach contact']
            # generate data frame from dictionary:
            df = generateDataframeCorrelationFromDic(data)

            #df = getDataFrameWT(df, 'C57BL/6J')

            plotCorrelationUsvEvents(ax=fig.add_subplot(gs[:, 15:20]),
                                     selectedEventList=selectedEventList,
                                     dataframe=df, cat='usvCorr', yMin=0, yMax=105, ageList=ageList,
                                     jitterValue=0.5,
                                     letter='D')

            ##################################################
            #plot the correlation of events with usvs
            selectedEventList = ['Oral-genital Contact', 'FollowZone Isolated', 'Train2']
            # open the json file to work on pre-computed data
            jsonFile = "dataCorrelationRefEvent.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file re-imported.")
            # generate data frame from dictionary:
            df = generateDataframeCorrelationFromDic(data)

            #df = getDataFrameKO(df, sex='female', age='3mo')
            print(df)

            plotCorrelationUsvEventsWithUSVWithKo(ax=fig.add_subplot(gs[:, 22:26]), selectedEventList=selectedEventList, dataframe=df,
                                     cat='usvCorr', yMin=0, yMax=105, ageList=ageList, jitterValue=0.3,
                                     letter='E')

            ################################
            ax = fig.add_subplot(gs[:, 27:30])

            ax.axis('off')
            malePatch = {}
            femalePatch = {}
            for age in ageList:
                malePatch[age] = matplotlib.lines.Line2D([0], [0], marker='v', color='w', markeredgecolor=getColorAge(age),
                                                           markerfacecolor=getColorAge(age), markeredgewidth=1.2, markersize=10,
                                                           label='B6 male {}'.format(age))
                femalePatch[age] = matplotlib.lines.Line2D([0], [0], marker='o', color='w', markeredgecolor=getColorAge(age), markerfacecolor=getColorAge(age), markeredgewidth=1.2,
                                                  markersize=10, label='B6 female {}'.format(age))

            koPatch = matplotlib.lines.Line2D([0], [0], marker='o', color='w', markeredgecolor=getColorKO(),
                                                       markerfacecolor=getColorKO(), markeredgewidth=1.2, markersize=10,
                                                       label='Shank3-/- female 3mo')

            handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'],
                       malePatch['7mo'], femalePatch['7mo'], koPatch]
            ax.legend(handles=handles, loc=(0.2, 0.1)).set_visible(True)

            fig.tight_layout()
            fig.savefig('figure_2_usv_usage.pdf', dpi = 300, bbox_inches="tight")
            fig.savefig('figure_2_usv_usage.jpg', dpi=300, bbox_inches="tight")
            plt.show()
            print("Job done.")

            break


        if answer == "3":
            yMinUsv = {'nbUsv': 0, 'nbBurst': 0, 'nbUsvBurst': 0, 'burstmeanDuration': 0, 'burstmeanInterval': 0}
            yMaxUsv = {'nbUsv': 24000, 'nbBurst': 5000, 'nbUsvBurst': 50, 'burstmeanDuration': 200, 'burstmeanInterval': 10}
            df = createDataframeFromJsonNumberUsvPerBurst(jsonFile='dataUsvDescription750.json', strainList=strainList, sexList=sexList, ageList=ageList)
            colSex = {'male': 'grey', 'female': 'black'}
            # Create general figure
            gs = gridspec.GridSpec(4, 3)
            # fig = plt.subplots(figsize=(24, 15), sharex=False, sharey=False)
            fig = plt.figure(figsize=(14, 17))

            # print(df)
            ###################################################
            # data wild-type
            dataframeWT = getDataFrameWT(df, strain='C57BL/6J')
            print('######')
            print(df['nbUsvBurst'])
            ax=fig.add_subplot(gs[0, 0:2])
            plotNumberUsvPerBurstWTBoxplot(ax=ax, dataframe=dataframeWT, yMinUsv=yMinUsv, yMaxUsv=yMaxUsv, letter='A', sexList=sexList, ageList=ageList)
            malePatch = {}
            femalePatch = {}
            for age in ageList:
                malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age),
                                                label='B6 male {}'.format(age))
                femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age),
                                                  label='B6 female {}'.format(age))

            handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'], malePatch['7mo'],
                       femalePatch['7mo']]
            ax.legend(handles=handles, loc=(0.23, 0.4)).set_visible(True)

            #########################################################
            # data KO number of USv per burst
            dataframeKO = getDataFrameKO(df, sex='female', age='3mo')
            ax=fig.add_subplot(gs[0, 2])
            plotNumberUsvPerBurstKOBoxplot(ax=ax, burstTrait='nbUsvBurst', dataframeKO=dataframeKO, yMinUsv=yMinUsv,
                                           yMaxUsv=yMaxUsv, letter='B')
            wtPatch = mpatches.Patch(edgecolor=colSex['female'], facecolor=getColorWT(), label='B6 female 3mo')
            koPatch = mpatches.Patch(edgecolor=colSex['female'], facecolor=getColorKO(), label='Shank3-/- female 3mo')
            handles = [wtPatch, koPatch]
            ax.legend(handles=handles, loc=(0.75, 0.7)).set_visible(True)

            # create model:
            print('nbUsvBurst')
            model = smf.mixedlm("nbUsvBurst ~ strain", dataframeKO, groups=dataframeKO['pair'])
            # run model:
            result = model.fit()
            # print summary
            print(result.summary())
            #plt.tight_layout()

            #number of USVs per burst accordind to the contexts for 3mo WT boxplot
            jsonFile = "burstTraitUsagePerEventContext.json"
            with open(jsonFile) as json_data:
                data = json.load(json_data)
            print("json file re-imported.")

            ax = fig.add_subplot(gs[1, 0:2])
            experimentsWT = getExperimentList(age='3mo', sex='female', genotype='WT')
            experimentsKO = getExperimentList(age='3mo', sex='female', strain='Shank3')
            plotBurstTraitsBoxplotsPerAgePerContextsPerGeno(letter='C', data=data, experimentsWT=experimentsWT,
                                                            experimentsKO=experimentsKO, ax=ax,
                                                            behavEventList=getFigureBehaviouralEvents(longList=False,
                                                                                                      withUrinate=False),
                                                            burstTrait='nbUSV', color=getColorKO())
            handles = [wtPatch, koPatch]
            ax.legend(handles=handles, loc=(0.25, 0.7)).set_visible(True)

            ####################################################
            #heatmap 3mo WT for nb of usv per burst:
            ax = fig.add_subplot(gs[1, 2])

            experiments = getExperimentList(age="3mo", sex="female", genotype="WT")
            #plotHeatmapBurstTraitUsagePerEventPerSet(experiments=experiments, genotype='WT', ax=ax, burstTrait='nbUSV', letter='D')
            
            '''
            ax.set_xticks( range( 0 , 9 ) )
            ax.set_yticks( range( 0 , 9 ) )
            ax.set_ylim( 0 , 9 )
            '''
            plotBurstTraitUsagePerEventContextPerSet2(experiments=experiments, genotype='WT', ax=ax, burstTrait='nbUSV', letter='D')
             

            '''# add scale on the heatmaps
            image = 'scale.jpg'
            imgPos = (3.5, -0.5)
            behavSchema = mpimg.imread(image)
            imgBox = OffsetImage(behavSchema, zoom=0.15)
            imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
            ax.add_artist(imageBox)'''
            


            # evolution of USV characteristics over the bursts
            letter = 'E'
            ax = fig.add_subplot(gs[2:4, 0:3])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.patch.set_facecolor('white')

            imgPos = (0.5, 0.5)
            image = 'figure6E.png'
            behavSchema = mpimg.imread(image)
            imgBox = OffsetImage(behavSchema, zoom=0.35)
            imageBox = AnnotationBbox(imgBox, imgPos, frameon=False)
            ax.add_artist(imageBox)

            ax.text(0, 1, letter, fontsize=20, horizontalalignment='center', color='black', weight='bold')

            plt.tight_layout()
            plt.show()
            fig.savefig('Figure_7_bursts_usage_750ms.pdf', dpi=300)
            print('Job done.')

            break

        if answer == "4":

            # Create general figure
            fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4.5))

            ####################################################
            # heatmap 5we WT for nb of usv per burst:
            ax = axes[0]
            experiments = getExperimentList(age="5we", sex="female", genotype="WT")
            plotBurstTraitUsagePerEventContextPerSet2(experiments=experiments, genotype='WT', ax=ax, burstTrait='nbUSV',
                                                      letter='A')
            ax.set_title('nb of USVs/burst', fontsize=16)
            ax.set_xlabel('B6 female 5we', fontsize=12, fontweight='bold', color=getColorAge('5we'))
            ax.set_facecolor('white')

            ####################################################
            # heatmap 7mo WT for nb of usv per burst:
            ax = axes[1]
            experiments = getExperimentList(age="7mo", sex="female", genotype="WT")
            plotBurstTraitUsagePerEventContextPerSet2(experiments=experiments, genotype='WT', ax=ax, burstTrait='nbUSV',
                                                      letter='B')
            ax.set_title('nb of USVs/burst', fontsize=16)
            ax.set_xlabel('B6 female 7mo', fontsize=12, fontweight='bold', color=getColorAge('7mo'))
            ax.set_facecolor('white')

            ####################################################
            # heatmap 3mo KO for nb of usv per burst:
            ax = axes[2]
            experiments = getExperimentList(age="3mo", sex="female", genotype="KO")
            plotBurstTraitUsagePerEventContextPerSet2(experiments=experiments, genotype='KO', ax=ax, burstTrait='nbUSV',
                                                      letter='C')

            ax.set_title('nb of USVs/burst', fontsize=16)
            ax.set_xlabel('Shank3-/- female 3mo', fontsize=12, fontweight='bold', color=getColorKO())
            ax.set_facecolor('white')


            plt.tight_layout()
            plt.show()
            fig.savefig('Suppl_Figure_heatmaps_bursts_750ms.pdf', dpi=300)
            fig.savefig('Suppl_Figure_heatmaps_bursts_750ms.jpg', dpi=300)
            print('Job done.')

            break


