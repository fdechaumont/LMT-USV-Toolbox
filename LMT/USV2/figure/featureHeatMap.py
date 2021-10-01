'''
Created on 1 juin 2021

@author: Fab

#inspired from:
#https://stackoverflow.com/questions/59381273/heatmap-with-circles-indicating-size-of-population

'''
from random import randint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

correspondanceList = [ (0.001, 1),
                      ( 0.01 , 0.9 ),
                      ( 0.05 , 0.8 ),
                      ( 0.1 , 0.2 )
                      ]
radiusScale = 3/7.0
effectSizeAmplitude = 1

def pValueToCircleSize( pValue ):
    
    size = 1
    for v in correspondanceList:
        if pValue >= v[0]:
            size = v[1]
    print( "pval to size: " , pValue, size )
    return size
    
def featureHeatMapPValLegend():
    fig, ax = plt.subplots(figsize=(3,6))
    y = 1
    M=1
    for v in correspondanceList:
        pValue = v[0]
        circle = plt.Circle((M , y), pValueToCircleSize(pValue) * radiusScale, color='grey', clip_on=False)
        plt.text(M, y + 0.6, "<=" + str(pValue), fontsize=10, ha="center", va="center")
        ax.add_patch(circle)
        y += 1.2
    ax.set_xlim(0,2)
    ax.set_ylim(0, 5)
    ax.invert_yaxis()
    plt.axis('off')
    plt.show()

def featureHeatMapEffectSizeLegend():
    fig, ax = plt.subplots(figsize=(3, 6))
    circles = []
    col = PatchCollection(circles, cmap="coolwarm")
    col.set_clim([-effectSizeAmplitude, effectSizeAmplitude])
    ax.add_collection(col)
    fig.colorbar(col)
    plt.show()

def featureHeatMap( dfEffectSize, dfPValue, ax , title=None,showLegend=False ):
    
    print("Data frame effect size:")            
    print( dfEffectSize )

    print("Data frame dfPValue:")            
    print( dfPValue )
    
    # labels
    
    xlabels = list( dfEffectSize.head() )
    ylabels = list( dfEffectSize.index )
        
    M= len( xlabels )
    N= len( ylabels )
    
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    print( x )
    print( y )
    

    circles = []
    colors = []
    s = []
    xx = 0
    for labX in xlabels:
        row = []
        yy=0
        for labY in ylabels:        
            effectSizeValue = dfEffectSize.loc[labY][labX]
            pValue = dfPValue.loc[labY][labX]
            row.append( dfEffectSize.loc[labY][labX] )
            print( xx , yy , effectSizeValue )
            circles.append ( Circle(( xx,yy ), radius=pValueToCircleSize(pValue) * radiusScale ) )
            colors.append( effectSizeValue )                    
            yy+=1
        xx+=1
        s.append(row)
    
    s = np.array( s )
    print( s )
    
    #fig, ax = plt.subplots( figsize=( len( xlabels ), len( ylabels ) ) )
    

    print ( circles )
    print( colors )
    col = PatchCollection(circles, array=np.array(colors), cmap="coolwarm" )
    # effect size scale
    col.set_clim([-effectSizeAmplitude, effectSizeAmplitude ])
    ax.add_collection(col)
        
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels , yticklabels= ""*len(ylabels) )
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which='minor')

    if title!=None:
        ax.set_title( title )

    ax.invert_yaxis()
    
    #plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    #plt.tight_layout()

    #FIXME
    #fig.colorbar(col)
    
    # draw p-value legend    
    if showLegend==True:
        y=0
        for v in correspondanceList:
            pValue= v[0]
            circle = plt.Circle( ( M+1, y ), pValueToCircleSize(pValue) * radiusScale , color='grey', clip_on=False)
            plt.text(M+1, y+0.6, "<="+str(pValue), fontsize=10, ha="center",va="center")
            ax.add_patch( circle )
            y+=1.2

    #plt.savefig( "test_"+ str( randint(1,1000))+".pdf" )
    return ax
    #plt.show()
    

    
    
def test():
    N = 10
    M = 11
    ylabels = ["".join(np.random.choice(list("PQRSTUVXYZ"), size=7)) for _ in range(N)]
    xlabels = ["".join(np.random.choice(list("ABCDE"), size=3)) for _ in range(M)]
    
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    s = np.random.randint(0, 180, size=(N,M))
    c = np.random.rand(N, M)-0.5
    
    fig, ax = plt.subplots()
    
    R = s/s.max()/2
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    print( circles )
    col = PatchCollection(circles, array=c.flatten(), cmap="RdYlGn")
    ax.add_collection(col)
    
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which='minor')
    
    fig.colorbar(col)
    plt.show()

if __name__ == '__main__':
    
    '''
    test()
    quit()
    '''
    
    colLabels = ["aaaaaaaaa","bbbbbbbbbb","cccccccccc","r","z","w","k"]
    rowLabels = ["eeeeeee","ffffff","gggggggggg","rze","er"]
    
    dataEffectSize= {}
    
    for col in colLabels:
        dataEffectSize[col] = []
        for animalType in rowLabels:
            dataEffectSize[col].append( randint( -20,20)/10 )
    
    dfEffectSize = pd.DataFrame( data=dataEffectSize, index = rowLabels )
    
    
    dataPValue= {}
    for event in colLabels:
        dataPValue[event] = []
        for animalType in rowLabels:
            dataPValue[event].append( randint(0,10)/100.0 )
    
    dfPValue = pd.DataFrame( data=dataPValue, index = rowLabels )
    
    
    featureHeatMap( dfEffectSize , dfPValue )
    
    
    
    
    
    
    
    
    