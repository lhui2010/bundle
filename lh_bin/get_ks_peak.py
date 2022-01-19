#!/usr/bin/env python

from sklearn.mixture import GaussianMixture
import logging
# import numpy as np
import pandas as pd
import sys
from iga.apps.base import sh

#Required python 3.8+ as default dict

class GausComp:
    def __init__(self, mean, sigma):
        self.mean  =str(mean)
        self.sigma =str(sigma)

class mGausComp:

    def __init__(self):
        self.dct = {}

    def append(self,mean,sigma):
        self.dct[mean] = GausComp(mean, sigma)

    def sort(self):
        new_dct = {}
        for k in sorted(self.dct):
            new_dct[k] = self.dct[k]
        self.dct = new_dct

    def get_result(self):
        self.sort()
        mean_res = ''
        sigma_res = ''
        for k in self.dct:
            mean_res += self.dct[k].mean + " "
            sigma_res += self.dct[k].sigma + " "
        return mean_res + "\t" + sigma_res


data=pd.read_table(sys.argv[1])
df=data.dropna(subset=['Ks'])
df=df[df.Ks<6]
GMM_input=df['Ks'].values
GMM_input=GMM_input.reshape(-1,1)
#res=GaussianMixture(n_components=3, covariance_type='spherical').fit(GMM_input)

if len(sys.argv) > 2:
    peak_num = int(sys.argv[2])
else:
    peak_num = 3

# Repeat 100 times
bic_dict = {}
with open(sys.argv[1] + ".KsGMM", 'w') as fh:
    for i in range(0,100):
        res=GaussianMixture(n_components=peak_num).fit(GMM_input)
        res_bic = res.bic(GMM_input)
        mean = []
        sigma = []
        result_comp = mGausComp()
        for m,s in zip(res.means_, res.covariances_):
            result_comp.append(float(m[0]), float(s))
        res_content = result_comp.get_result()
        result_line = "{}\t{}".format(res_bic, res_content)
        fh.write(result_line + "\n")
        bic_dict[res_bic] = res_content

#print best with lowest bic
#format: means seperated with ' ' \t sigma seperated with ' '
print(bic_dict[sorted(bic_dict)[0]])








##Deprecated
# mylist = []
# with open(sys.argv[1]) as fh:
#     for line in fh:
#         mylist.append(line.rstrip())
# mylist = np.array(mylist)
# mylist = mylist.reshape(-1, 1)
# for s in res.covariances_:
#     sigma.append(s)
# print(sigma)
# res=GaussianMixture(n_components=3, covariance_type='spherical').fit(GMM_input)
# mean
# print(res.means_.tolist())
# print(res.covariances_.tolist())

