'''
Created on 1 oct. 2021
Modified by Elodie on 23 Jan 2026

@author: Fab
'''


def fileBiggerSplit( fileName ):

    fileName = fileName.replace("-","_")
    fileName = fileName.replace(".","_")    
    result =""
    for s in fileName.split("_"):
        if len( s ) > len( result ):
            result = s
    return result


def getDataFileMatch( dataFiles, number, extension = ".txt" ):
    for file in dataFiles:
        if not extension in file:
            continue
        if number in file:
            return file
    return None


if __name__ == '__main__':
    pass