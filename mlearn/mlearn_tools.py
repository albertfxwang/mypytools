#!/usr/bin/env python

import numpy as np
import sklearn
from sklearn import cross_validation

def cross_val_score(model, trainMatrix, trainLabels, cv=5):
   scores = cross_validation.cross_val_score(model, trainMatrix, trainLabels, 
                                             cv=cv)
   print "Accuracy: %0.2f +/- %0.2f" % (scores.mean(), scores.std()*2)

def predict_labels(model, testMatrix, testNames, output=""):
   testLabels = model.predict(testMatrix)
   if type(testLabels[0]) == type('0'):
      testLabels = [chr(x) for x in testLabels]
   testNames = [int(os.path.splitext(x)[0]) for x in testNames]
   if len(output):
      with open(output, 'wb') as f:
         f.write("ID,Class\n")
         for i in range(len(testNames)):
            f.write("%d,%s\n" % (testNames[i], testLabels[i]))
   print "Done."