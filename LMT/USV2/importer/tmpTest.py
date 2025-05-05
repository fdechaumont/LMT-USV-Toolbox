'''
Created on 21 mars 2022

@author: Fab
'''
from LMT.USV2.importer.Voc import Voc
import pickle


if __name__ == '__main__':
    print("testing")
    voc = Voc()
    print( voc )
    rf = pickle.load( open("trainingSet.bin", 'rb') )
    