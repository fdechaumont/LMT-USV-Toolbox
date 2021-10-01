'''
Created on 20 janv. 2020

@author: Elodie
'''


from scipy.stats import kruskal, ttest_ind

from lmtanalysis.FileUtil import extractPValueFromLMMResult
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf

from LMT.USV2.figure.figUtil import *
from LMT.USV2.figure.figParameter import *
from LMT.USV2.figure.vocUtil import *
from LMT.USV2.experimentList.experimentList import getAllExperimentList
import sqlite3

from lmtanalysis.Event import EventTimeLine
import json
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from lmtanalysis.Util import getMinTMaxTInput
import pandas
import string
import matplotlib.patches as mpatches
from LMT.USV2.figure.figUtil import sexList, strainList, yMinVar, yMaxVar
from LMT.USV2.figure.figParameter import getFigureVocTraits
from lmtanalysis.Animal import AnimalPool
from LMT.USV2.lib.USVUtil import getStrainAgeSexPairGenoPerFile


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


def plotAcousticVarWTBoxplot(ax, dataframe, yMinUsv, yMaxUsv, variable, letter):
    xPos = {'5we': 0, '3mo': 1, '7mo': 2}
    yLabel = getFigureLabelTrait(variable)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 9)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='age', y=variable, hue='sex', hue_order=['male', 'female'], data=dataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    ax.set_xticklabels(['5we', '3mo', '7mo'], rotation=0, fontsize=16, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=18)
    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)
    bp.legend().set_visible(False)

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

    '''
                    data=dfWT.loc[:, ['age', 'sex', 'pair', str(var)]]
                    formula = 'var ~ C(age) + C(sex) + C(age):C(sex)'
                    model = ols(formula, data).fit()
                    aov_table = anova_lm(model, typ=2)
                    '''
    n = 0
    for ageClass in ageList:
        valM = dataframe[variable][(dataframe['age'] == ageClass) & (dataframe['sex'] == 'male')]
        valF = dataframe[variable][(dataframe['age'] == ageClass) & (dataframe['sex'] == 'female')]
        nMales = len(valM)
        nFemales = len(valF)
        meanM = np.mean(valM)
        meanF = np.mean(valF)
        print('Tests between sexes: females {} voc, males {} voc'.format(nFemales, nMales))

        # Mixed model: variable to explain: value; fixed factor = genotype; random effect: group
        dfTestAllAges = dataframe[['age', 'sex', 'pair', 'strain', variable]]
        dfTest = dfTestAllAges[dfTestAllAges['age'] == ageClass]
        dfTest.rename(columns={'age': 'age', 'sex': 'sex', 'pair': 'pair', variable: 'value'}, inplace=True)
        # create model:
        model = smf.mixedlm("value ~ sex", dfTest, groups=dfTest['pair'])
        # run model:
        result = model.fit()
        # print summary
        print(variable, ageClass)
        print(result.summary())
        p, sign = extractPValueFromLMMResult(result=result, keyword='male')

        '''U, p = ttest_ind(
            dataframe[variable][(dataframe['age'] == ageClass) & (dataframe['sex'] == 'male')],
            dataframe[variable][(dataframe['age'] == ageClass) & (dataframe['sex'] == 'female')])
        print(variable, ' ', ageClass, ' males ', meanM, ' females ', meanF, ' U = ', U, 'p = ', p,
              getStarsFromPvalues(p, 3 * len(acousticVariables)))'''

        # add p-values on the plot
        correction = 3 * len(acousticVariables)
        if p * correction >= 0.05:
            stars = getStarsFromPvalues(p, 1)
        elif p * correction < 0.05:
            stars = getStarsFromPvalues(p, 1) + '°'
        ax.text(xPos[ageClass],
                yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]), stars,
                fontsize=20, horizontalalignment='center', color='black', weight='bold')
        n += 1
    return bp


def plotAcousticVarKOBoxplot(ax, dataframe, yMinUsv, yMaxUsv, variable, letter):
    yLabel = getFigureLabelTrait(variable)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2, 4, 5, 7, 8]
    ax.set_xticks(xIndex)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 9)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-0.8, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='black', markeredgecolor='black')
    bp = sns.boxplot(x='age', y=variable, hue='strain', data=dataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    ax.set_xticklabels(['3mo'], rotation=0, fontsize=16, horizontalalignment='center')
    ax.set_ylabel(yLabel, fontsize=18)
    ax.set_xlabel(ageList, fontsize=18)

    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)
    bp.set_ylabel(yLabel, fontsize=15)
    bp.legend().set_visible(False)

    colorList = [getColorWT(), getColorKO()]
    edgeList = ['black', 'black']
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))


    valB6 = dataframe[variable][(dataframe['strain'] == 'C57BL/6J') & (dataframe['age'] == '3mo') & (dataframe['sex'] == 'female')]
    valShank3 = dataframe[variable][(dataframe['strain'] == 'Shank3') & (dataframe['age'] == '3mo') & (dataframe['sex'] == 'female')]
    nB6 = len(valB6)
    nShank3 = len(valShank3)
    meanB6 = np.mean(valB6)
    meanShank3 = np.mean(valShank3)
    print('Tests between B6 vs Shank3: B6 {} voc, Shank3 {} voc'.format(nB6, nShank3))

    # Mixed model: variable to explain: value; fixed factor = genotype; random effect: pair
    dfTest = dataframe[['age', 'sex', 'pair', 'strain', variable]]
    dfTest.rename(columns={'age': 'age', 'sex': 'sex', 'pair': 'pair', variable: 'value'}, inplace=True)
    # create model:
    model = smf.mixedlm("value ~ strain", dfTest, groups='pair')
    # run model:
    result = model.fit()
    # print summary
    print('################')
    print(variable)
    print(result.summary())
    p, sign = extractPValueFromLMMResult(result=result, keyword='Shank3')

    '''T, p = ttest_ind(
        dataframe[variable][(dataframe['strain'] == 'C57BL/6J') & (dataframe['age'] == '3mo') & (dataframe['sex'] == 'female')],
        dataframe[variable][(dataframe['strain'] == 'Shank3') & (dataframe['age'] == '3mo') & (dataframe['sex'] == 'female')])
    print(variable, ' B6 ', meanB6, ' vs Shank3 ', meanShank3, ': U = ', T, 'p = ', p,
          getStarsFromPvalues(p, 1 * len(acousticVariables)))'''
    # add p-values on the plot
    correction = 1 * len(acousticVariables)
    if p * correction >= 0.05:
        stars = getStarsFromPvalues(p, 1)
    elif p * correction < 0.05:
        stars = getStarsFromPvalues(p, 1) + '°'
    ax.text(0,
            yMaxVar[variable] - 0.06 * (yMaxVar[variable] - yMinVar[variable]), stars,
            fontsize=16, horizontalalignment='center', color='black', weight='bold')

    return bp


def plotAcousticVarBoxplotDiffGeno(ax, dataframe, yMinUsv, yMaxUsv, variable, genoList, letter):
    yLabel = getFigureLabelTrait(variable)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, 2)
    ax.set_ylim(yMinUsv[variable], yMaxUsv[variable])
    ax.tick_params(axis='y', labelsize=14)

    ax.text(-1.5, yMaxUsv[variable] + 0.06 * (yMaxUsv[variable] - yMinUsv[variable]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    meanprops = dict(marker='o', markerfacecolor='white', markeredgecolor='black')
    bp = sns.boxplot(x='geno', y=variable, data=dataframe, ax=ax, showfliers = False, showmeans=True, notch=True, meanprops=meanprops, width=0.4, dodge=True)
    ax.set_xticklabels(genoList, rotation=45, fontsize=12, horizontalalignment='right')

    bp.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0)

    bp.set(xlabel=None)

    bp.set_ylabel(yLabel, fontsize=15)
    bp.legend().set_visible(False)

    colorList = [getColorGeno(genoList[0]), getColorGeno(genoList[1])]
    edgeList = ['black', 'black']
    n = 0
    for box in bp.artists:
        box.set_facecolor(colorList[n])
        box.set_edgecolor(edgeList[n])
        n += 1
    # Add transparency to colors
    for box in bp.artists:
        r, g, b, a = box.get_facecolor()
        box.set_facecolor((r, g, b, .7))

    # Mixed model: variable to explain: value; fixed factor = genotype; random effect: group
    dataframeToTest = dataframe.loc[:,['strain', 'sex', 'geno', 'pair', variable]]
    dataframeToTest.rename(columns={variable: 'value'}, inplace=True)
    # create model:
    model = smf.mixedlm("value ~ geno", dataframeToTest, groups=dataframeToTest['pair'])
    # run model:
    result = model.fit()
    # print summary
    print(variable)
    print(result.summary())
    p, sign = extractPValueFromLMMResult(result=result, keyword='WT')
    ax.text(0.5, yMaxUsv[variable] - 0.06 * (yMaxUsv[variable] - yMinUsv[variable]),
            getStarsFromPvalues(pvalue=p, numberOfTests=1), fontsize=20)


def plotAcousticVarAgeEffect(ax, dfWT, var, ageList, yMinVarAge, yMaxVarAge, colSex, letter):
    meanVal = {}
    semVal = {}
    for sex in ['male', 'female']:
        meanVal[sex] = {}
        semVal[sex] = {}
        for ageClass in ageList:
            val = dfWT[var][(dfWT['age'] == ageClass) & (dfWT['sex'] == sex)]
            n = len(val)
            meanVal[sex][ageClass] = np.mean(list(val))
            # semVal[sex][ageClass] = np.std(list(val))/np.sqrt(n) #sem
            semVal[sex][ageClass] = np.std(list(val))  # std

    yLabel = getFigureLabelTrait(var)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xIndex = [1, 2, 3]
    ax.set_xticks(xIndex)

    ax.set_ylabel(yLabel, fontsize=16)
    ax.xaxis.set_tick_params(direction="in")
    ax.yaxis.set_tick_params(direction="in")
    ax.set_xlim(0, len(ageList) + 1)
    ax.set_ylim(yMinVarAge[var], yMaxVarAge[var])
    ax.tick_params(axis='y', labelsize=14)

    xF = np.array(xIndex) + 0.1
    xM = np.array(xIndex) - 0.1

    ax.text(-0.5, yMaxVarAge[var] + 0.06 * (yMaxVarAge[var] - yMinVarAge[var]), letter,
            fontsize=20, horizontalalignment='center', color='black', weight='bold')
    ax.errorbar(xF, [meanVal['female']['5we'], meanVal['female']['3mo'], meanVal['female']['7mo']],
                yerr=[semVal['female']['5we'], semVal['female']['3mo'], semVal['female']['7mo']], marker='o', capsize=4,
                color=colSex['female'], ecolor=colSex['female'])
    ax.errorbar(xM, [meanVal['male']['5we'], meanVal['male']['3mo'], meanVal['male']['7mo']],
                yerr=[semVal['male']['5we'], semVal['male']['3mo'], semVal['male']['7mo']], marker='v', capsize=4,
                color=colSex['male'], ecolor=colSex['male'])
    ax.set_xticklabels(ageList, rotation=0, fontsize=16,
                       horizontalalignment='center')
    # ax.plot(xIndex, [meanVal['female']['5we'], meanVal['female']['3mo'], meanVal['female']['7mo']], color='red', zorder=2)
    # ax.plot(xIndex, [meanVal['male']['5we'], meanVal['male']['3mo'], meanVal['male']['7mo']], color='blue', zorder=2)
    # ax.scatter(xIndex, [meanVal['female']['5we'], meanVal['female']['3mo'], meanVal['female']['7mo']], c='red', zorder=1)
    # ax.scatter(xIndex, [meanVal['male']['5we'], meanVal['male']['3mo'], meanVal['male']['7mo']], c='blue', zorder=1)

    posSex = {'male': [1.75, 2.75], 'female': [1.25, 2.25]}
    for sex in sexList:
        K, p = kruskal(
            dfWT[var][(dfWT['age'] == '5we') & (dfWT['sex'] == sex)],
            dfWT[var][(dfWT['age'] == '3mo') & (dfWT['sex'] == sex)],
            dfWT[var][(dfWT['age'] == '7mo') & (dfWT['sex'] == sex)])
        print(var, ' ', sex, ' K = ', K, 'p = ', p, getStarsFromPvalues(p, 1 * len(acousticVariables)))
        if p < 0.05:
            U1, p1 = ttest_ind(
                dfWT[var][(dfWT['age'] == '5we') & (dfWT['sex'] == sex)],
                dfWT[var][(dfWT['age'] == '3mo') & (dfWT['sex'] == sex)])
            print(var, sex, ' 5we versus 3mo: U = ', U1, 'p = ', p1,
                  getStarsFromPvalues(p1, 2 * len(acousticVariables)))
            U2, p2 = ttest_ind(
                dfWT[var][(dfWT['age'] == '3mo') & (dfWT['sex'] == sex)],
                dfWT[var][(dfWT['age'] == '7mo') & (dfWT['sex'] == sex)])
            print(var, sex, ' 3mo versus 7mo: U = ', U2, 'p = ', p2,
                  getStarsFromPvalues(p1, 2 * len(acousticVariables)))
            correction = 1 * len(acousticVariables)
            if p1 * correction >= 0.05:
                stars1 = getStarsFromPvalues(p1, 1)
            elif p1 * correction < 0.05:
                stars1 = getStarsFromPvalues(p1, 1) + '°'
            ax.text(posSex[sex][0],
                    yMaxVarAge[var] - 0.2 * (yMaxVarAge[var] - yMinVarAge[var]),
                    stars1,
                    fontsize=16, horizontalalignment='center', color=colSex[sex])
            if p2 * correction >= 0.05:
                stars2 = getStarsFromPvalues(p2, 1)
            elif p2 * correction < 0.05:
                stars2 = getStarsFromPvalues(p2, 1) + '°'
            ax.text(posSex[sex][1],
                    yMaxVarAge[var] - 0.2 * (yMaxVarAge[var] - yMinVarAge[var]),
                    stars2,
                    fontsize=16, horizontalalignment='center', color=colSex[sex])


if __name__ == '__main__':
    '''
    This codes extracts the acoustic variables of USV according to the contexts of emission. Data are stored in a json file for further computation.
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
    print(colNames)

    #yMinVar, yMaxVar, yMinVarAge and yMaxVarAge imported from vocUtil


    while True:

        question = "Do you want to:"
        question += "\n\t [e]xtract acoustic variables on the whole set of data (save json file)?"
        question += "\n\t run the [pca] analysis on the whole set of C57BL/6J data (from stored json file)?"
        question += "\n\t [pv] plot and analyse raw acoustic variables in C57BL/6J mice (from stored json file)?"
        question += "\n\t [page] simple plot and analyse raw acoustic variables for age in C57BL/6J mice (from stored json file)?"
        question += "\n\t [pko] plot and analyse raw acoustic variables in C57BL/6J and Shank3 KO mice (from stored json file)?"
        question += "\n"
        answer = input(question)

        if answer == "e":
            experiments = getAllExperimentList()

            tmin, tmax = getMinTMaxTInput()

            # initialise dictionaries to store data
            data = {}
            for strain in strainList:
                data[strain] = {}
                for sex in sexList:
                    data[strain][sex] = {}
                    for age in ageList:
                        data[strain][sex][age] = {}

            dataDf = {}
            for col in colNames:
                dataDf[col] = []

            #extract data from files
            for exp in experiments:
                print('exp: ', exp.getFullName())

                file = exp.file
                print( file )
                connection = sqlite3.connect( file )
                pool = AnimalPool()
                pool.loadAnimals(connection)
                strainAgeSexPairGenoFile = getStrainAgeSexPairGenoPerFile( connection )

                strainFile = strainAgeSexPairGenoFile[0]
                ageFile = strainAgeSexPairGenoFile[1]
                sexFile = strainAgeSexPairGenoFile[2]
                pairFile = strainAgeSexPairGenoFile[3]
                genoFile = strainAgeSexPairGenoFile[4]

                #initialise the datalist for the specific file:
                data[strainFile][sexFile][ageFile][pairFile] = {}
                for variable in acousticVariables:
                    data[strainFile][sexFile][ageFile][pairFile][variable] = []

                #clean the voc timeline using the "excluded" metadata
                print("Cleaning voc with excluded metadata")
                vocTimeLine = EventTimeLine( connection, "Voc", minFrame = tmin, maxFrame = tmax, loadEventIndependently = True )
                
                print("remaining voc events: ", len(vocTimeLine.eventList))

                #Select the vocalisations:
                for vocEvent in vocTimeLine.getEventList():
                    if vocEvent.metadata['excluded'] == True:
                        print('excluded!')
                        continue
                    else:
                        dataDf['strain'].append(strainFile)
                        dataDf['sex'].append(sexFile)
                        dataDf['age'].append(ageFile)
                        dataDf['pair'].append(pairFile)
                        for variable in acousticVariables:
                            data[strainFile][sexFile][ageFile][pairFile][variable].append( vocEvent.metadata[variable] )
                            dataDf[variable].append(vocEvent.metadata[variable])

                connection.close()

            #Create a json file to store the computation
            with open("dataAcousticAnalysisAllUsvs.json", 'w') as fp:
                json.dump( dataDf, fp, indent=4 )
            print("json file with acoustic measurements created.")

            break


        if answer == 'pca':

            #create a panda dataframe with acoustic measures for all voc
            with open("dataAcousticAnalysisAllUsvs.json") as json_data:
                dataDf = json.load(json_data)

            print("json file for acoustic variables re-imported.")
            df = pd.DataFrame( dataDf, columns= colNames )
            print(df.head( n = 10 ))
            print( 'Dataframe created.')

            #select the data only for C57BL/6J mice
            dfWT = df.loc[df['strain']=='C57BL/6J', acousticVariables]

            #prepare ACP: verify version of scikit-learn:
            print('version sklearn: ', sklearn.__version__)

            # instanciation
            sc = StandardScaler()

            # transformation – center-reduce
            dfWTTrans = sc.fit_transform(dfWT)
            #print(dfWTTrans)

            # mean and std verification
            print('means: ', np.mean(dfWTTrans,axis=0))
            print( 'std: ', np.std(dfWTTrans, axis=0, ddof=0))

            # instanciation
            acp = PCA(svd_solver='full')
            # affichage des paramètres
            print(acp)

            # calculs
            coord = acp.fit_transform(dfWTTrans)
            #nombre de composantes calculées
            print(acp.n_components_) # nb of variables
            n = acp.n_components_
            #explained variance
            print(acp.explained_variance_)
            #corrected value
            eigval = (n-1)/n*acp.explained_variance_
            print('corrected values for explained variance: ', eigval)

            # proportion of explained variance
            print('proportion of explained variance: ', acp.explained_variance_ratio_)

            # cumul of explained variance
            plt.plot(np.arange(1,n+1,1),np.cumsum(acp.explained_variance_ratio_))
            plt.title("Explained variance vs. # of factors")
            plt.ylabel("Cumsum explained variance ratio")
            plt.xlabel("Factor number")
            plt.show()

            #determine the contribution of each variable to the components
            sqrt_eigval = np.sqrt(eigval)
            corvar = np.zeros((n, n))
            for k in range(n):
                corvar[:, k] = acp.components_[k, :] * sqrt_eigval[k]

            # afficher la matrice des corrélations variables x facteurs
            print(corvar)
            # on affiche pour les 6 premiers axes
            contribVar = pandas.DataFrame({'id':acousticVariables,'COR_1':corvar[:,0],'COR_2':corvar[:,1],'COR_3':corvar[:,2],'COR_4':corvar[:,3],'COR_5':corvar[:,4],'COR_6':corvar[:,5]})
            print(contribVar)

            contribVar.to_excel(r'C:\Users\eye\Documents\GitHub\USV-analysis\LMT\USV\testVocTraitUsagePerEventType\contribVarPCA.xlsx', index=False, header=True)

            print("Job done.")
            break

        if answer == 'pv':
            # create a panda dataframe with acoustic measures for all voc
            with open("dataAcousticAnalysisAllUsvs.json") as json_data:
                dataDf = json.load(json_data)

            print("json file for acoustic variables re-imported.")
            df = pd.DataFrame(dataDf, columns=colNames)
            print(df.head(n=10))
            print('Dataframe created.')
            #xPos = {'5we': 1.5, '3mo': 4.5, '7mo': 7.5}
            xPos = {'5we': 0, '3mo': 1, '7mo': 2}
            # select the data only for C57BL/6J mice
            dfWT = df.loc[df['strain'] == 'C57BL/6J', :]
            row = 0
            col = 0
            letterList = list(string.ascii_uppercase)
            k = 0
            fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(24, 12))

            for var in acousticVariables:
                print('################################')
                print(var)
                ax = axes[row][col]
                plotAcousticVarWTBoxplot( ax=ax, dataframe=dfWT, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[k])


                if col < 3:
                    col += 1
                else:
                    col = 0
                    row += 1
                k += 1

            ##############################
            ax = axes[0, 0]

            malePatch = {}
            femalePatch = {}
            for age in ageList:
                malePatch[age] = mpatches.Patch(edgecolor='grey', facecolor=getColorAge(age),
                                                label='male {}'.format(age))
                femalePatch[age] = mpatches.Patch(edgecolor='black', facecolor=getColorAge(age),
                                                  label='female {}'.format(age))

            handles = [malePatch['5we'], femalePatch['5we'], malePatch['3mo'], femalePatch['3mo'], malePatch['7mo'],
                       femalePatch['7mo']]
            ax.legend(handles=handles, loc=(0.05, 0.36)).set_visible(True)

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_var_raw_male_female_all_ages.pdf', dpi=300)
            fig.savefig('fig_var_raw_male_female_all_ages.jpg', dpi=300)
            print('Job done.')
            break


        if answer == 'page':
            # create a panda dataframe with acoustic measures for all voc
            with open("dataAcousticAnalysisAllUsvs.json") as json_data:
                dataDf = json.load(json_data)

            print("json file for acoustic variables re-imported.")
            df = pd.DataFrame(dataDf, columns=colNames)
            print(df.head(n=10))
            print('Dataframe created.')
            xPos = {'5we': 1, '3mo': 2, '7mo': 3}
            # select the data only for C57BL/6J mice
            dfWT = df.loc[df['strain'] == 'C57BL/6J', :]
            row = 0
            col = 0
            letterList = list(string.ascii_uppercase)
            colSex = {'male': 'grey', 'female': 'black'}
            k = 0
            fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(24, 12))

            for var in acousticVariables:
                print('################################')
                print(var)
                ax = axes[row][col]
                plotAcousticVarAgeEffect( ax=ax, dfWT=dfWT, var=var, yMinVarAge=yMinVarAge, yMaxVarAge=yMaxVarAge, ageList=ageList, colSex=colSex, letter=letterList[k] )

                if col < 3:
                    col += 1
                else:
                    col = 0
                    row += 1
                k += 1

            ax = axes[0, 0]
            malePatch = mpatches.Patch(edgecolor=colSex['male'], facecolor='white', label='male')
            femalePatch = mpatches.Patch(edgecolor=colSex['female'], facecolor='white', label='female')
            handles = [malePatch, femalePatch]
            ax.legend(handles=handles, loc=(0.05, 0.6)).set_visible(True)

            plt.tight_layout()
            plt.show()
            fig.savefig('fig_var_age effect_male_female.pdf', dpi=300)
            fig.savefig('fig_var_age effect_male_female.jpg', dpi=300)
            print('Job done.')
            break


        if answer == 'pko':
            # create a panda dataframe with acoustic measures for all voc
            with open("dataAcousticAnalysisAllUsvs.json") as json_data:
                dataDf = json.load(json_data)

            print("json file for acoustic variables re-imported.")
            df = pd.DataFrame(dataDf, columns=colNames)
            print(df.head(n=10))
            print('Dataframe created.')
            # select the data only for C57BL/6J and Shank3 3mo female mice
            dfKO = df.loc[(df['sex'] == 'female') & (df['age'] == '3mo'), :]
            row = 0
            col = 0
            letterList = list(string.ascii_uppercase)
            k = 0
            fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(18, 12))

            for var in acousticVariables:
                print(var)
                ax = axes[row][col]
                plotAcousticVarKOBoxplot( ax=ax, dataframe=dfKO, yMinUsv=yMinVar, yMaxUsv=yMaxVar, variable=var, letter=letterList[k])

                if col < 3:
                    col += 1
                else:
                    col = 0
                    row += 1
                k += 1

            ax = axes[0,0]
            wtPatch = mpatches.Patch(edgecolor='black', facecolor=getColorWT(), label='B6 female')
            koPatch = mpatches.Patch(edgecolor='black', facecolor=getColorKO(), label='Shank3 female')
            handles = [wtPatch, koPatch]
            ax.legend(handles=handles, loc=(0.05, 0.7)).set_visible(True)

            plt.tight_layout()
            plt.show()

            fig.savefig('fig_var_raw_female_B6_KO.pdf', dpi=300)
            fig.savefig('fig_var_raw_female_B6_KO.jpg', dpi=300)
            print('Job done.')
            break