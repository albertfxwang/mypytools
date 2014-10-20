import numpy as np

class Gini(object):
   # calculate the Gini coefficient, following different estimators
   def __init__(self, data):
      self.data = np.array(data).astype('float')
      self.data_sorted = np.sort(data)  # sort data in increasing order

   def rel_mean_diff(self):
      # follow the formula here:
      # http://mathworld.wolfram.com/GiniCoefficient.html
      # Can use unordered data.
      # first, calculate the sum of every possible absolute pair difference
      d = np.roll(self.data, 1)  # shift all element by 1 to the right
      pairsum = 0.
      for i in range(len(d) - 1):
         p = np.abs(self.data - d)
         pairsum += p.sum()
         d = np.roll(d, 1)
      ndata = len(self.data)
      g = pairsum / (2. * ndata**2 * np.mean(self.data))
      return g

   def formula1(self):
      # requires that data be sorted in non-decreasing order
      # the first formula on the Wikipedia page:
      # http://en.wikipedia.org/wiki/Gini_coefficient
      ndata = len(self.data_sorted)
      q = np.sum((ndata + 1 - np.arange(1, ndata+1)) * self.data_sorted)
      q = q / np.sum(self.data_sorted)
      g = 1. / ndata * (ndata + 1. - 2. * q)
      return g