import numpy as np
from Bio import SeqIO
import re
import sys
import glob
import statistics
import pywt
#from sympy import fft es muy lento
from scipy.fftpack import fft
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn import svm
from sklearn.model_selection import train_test_split
import joblib
import os

from train import descriptor

