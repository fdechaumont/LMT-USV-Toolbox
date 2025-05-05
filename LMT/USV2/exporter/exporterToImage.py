'''
Created on 23 janv. 2025

@author: Fab
'''
import glob
import os
from LMT.USV2.exporter.MachineLearningTester import MachineLearningTester
from LMT.USV2.importer.Wav import Wav




    



if __name__ == '__main__':
    
    
    path = "D:/collab voc data/converted/render/"
    
    m = MachineLearningTester()
    
    files = glob.glob("D:/collab voc data/converted/*.wav")
    
    for file in files:
        
        print ( "current file: " , files.index( file ) , " / " , len ( files ))
        
        vocList = m.process( file )
        print( len( vocList ) )
        if len( vocList ) > 4:
            
            w = Wav( file )
            fileName= os.path.basename( file )
            fileName = path+fileName+".jpg"
            print( fileName )
            try:            
                w.saveSpectrumImage( fileName , minT=None, maxT = None , text= None , textColor= "black", save=True , vmin= -5 , vmax = 10 )
            except:
                pass
            
            
            
            
            
            
        
    
     
    