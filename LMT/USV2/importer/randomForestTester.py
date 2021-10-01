'''
Created on 24 sept. 2018

@author: Fab

'''
import unittest

# Load the library with the iris dataset
from sklearn.datasets import load_iris

# Load scikit's random forest classifier library
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree

# Load pandas
import pandas as pd

# Load numpy
import numpy as np
from sklearn.model_selection._validation import cross_val_score

class RandomForestUnit:
    
    def __init__( self , df , testCLF=True , showAsATree = False ):
        
        # Create a new column that for each row, generates a random number between 0 and 1, and
        # if that value is less than or equal to .75, then sets the value of that cell as True
        # and false otherwise. This is a quick and dirty way of randomly assigning some rows to
        # be used as the training data and some as the test data.
        df['is_train'] = np.random.uniform(0, 1, len(df)) <= .4

        
        # Create two new dataframes, one with the training rows, one with the test rows
        train, test = df[df['is_train']==True], df[df['is_train']==False]
        
        # Show the number of observations for the test and training dataframes
        #print('Number of observations in the training data:', len(train))
        #print('Number of observations in the test data:',len(test))
        
        # Create a list of the feature column's names
        # assumes target is the last column
        features = df.columns[:-2]
        
        # View features
        print ("Features:")
        print ( features )
        
        # train['target'] contains the actual target names. Before we can use it,
        # we need to convert each species name into a digit. So, in this case there
        # are three species, which have been coded as 0, 1, or 2.
        y , self.targetCorrespondanceNames = pd.factorize(train['target'])
        # View target
        #print ("indexes" , y, "target indexes:" , self.targetCorrespondanceNames )
        
        # Create a random forest Classifier. By convention, clf means 'Classifier'
        if showAsATree == True :
            self.clf = tree.DecisionTreeClassifier()
        else:
            self.clf = RandomForestClassifier(n_jobs=2, random_state=100, n_estimators=1000 )
        
        
        # Train the Classifier to take the training features and learn how they relate
        # to the training y (the targets)
        #print ("Train features")
        #print ( train[features] )
        
        self.clf.fit(train[features], y)
        
        '''
        if showAsATree == True :
            import graphviz             
            dot_data = tree.export_graphviz(self.clf, out_file=None , feature_names= features, class_names=self.targetCorrespondanceNames, filled=True, rounded=True, special_characters=True  ) 
            graph = graphviz.Source(dot_data) 
            graph.render("Decision Tree") 
        '''
        '''
        dot_data = tree.export_graphviz(clf, out_file=None, 
                         feature_names=iris.feature_names,  
                         class_names=iris.target_names,  
                         filled=True, rounded=True,  
                         special_characters=True) 
        '''
        
        if testCLF == False:
             return

        # Apply the Classifier we trained to the test data (which, remember, it has never seen before)
        self.clf.predict(test[features])
        
        # View the predicted probabilities of the first 10 observations
        #print ( clf.predict_proba(test[features])[0:10] )
        
        # Create actual names for the target
        preds = self.clf.predict(test[features])
        
        # View the PREDICTED species for the first five observations
        #print( "predicted (head):")
        #print ( preds[0:5] )
        
        # View the ACTUAL targets for the first five observations
        #print( test['target'].head() )
        
        # Create confusion matrix
        print( "Confusion matrix:")
        
        #confusion_matrix(y_true, y_pred)
        
        print( pd.crosstab(test['target'], preds, rownames=['Actual targets'], colnames=['Predicted targets']) )
        
        # View a list of the features and their importance scores
        # print ( list(zip(train[features], clf.feature_importances_)) )
    
        # another view:
        
        # Get numerical feature importances
        importances = list(self.clf.feature_importances_)

        # List of tuples with variable and importance
        feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(features, importances)]
        
        # Sort the feature importances by most important first
        feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
        
        # Print out the feature and importances 
        [print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances]

        print( "computing cross validation score")
        # the std of accuracy could be reduced by launching more tests. (to test)
        scores = cross_val_score(self.clf, train[features], y, cv=10, n_jobs=-1 )
        print ( scores )
        print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() ))
        
        self.accuracy = scores.mean()
        self.accuracyError = scores.std()

    def getTargetCorrespondanceNames(self):
        return self.targetCorrespondanceNames

class TestRandomForestUnit ( unittest.TestCase ):
    
    def test_Species(self):

        iris = load_iris()

        pd.set_option('display.expand_frame_repr', False)

        # Create a dataframe with the four feature variables
        df = pd.DataFrame(iris.data, columns=iris.feature_names)
        
        df['target'] = pd.Categorical.from_codes(iris.target, iris.target_names)
        
        print( df )

        rf = RandomForestUnit( df )

        self.assertEqual( True, True )
    
    