#!/usr/bin/env python

import numpy as np

"""
Some simple statistical functions.
"""

def consistent_sigma(val1, sig1, val2, sig2):
   """
   Given the values and erros for both values, test if they are 
   consistent within the errors.
   """
   assert (sig1 > 0) & (sig2 > 0)
   if val1 >= val2:
      return (val2 + sig2) >= (val1 - sig1)
   else:
      return (val1 + sig1) >= (val2 - sig2)
      