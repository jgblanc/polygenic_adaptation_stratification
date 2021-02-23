# This script generates Test vector and calculates Tm 

# Import modules
import argparse
from pandas_plink import read_plink
import pandas as pd
import dask
import dask_ml.preprocessing
import numpy as np

# Parse Inputs 
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")
req_grp.add_argument("--inpre","-i",dest="inpre",help="input plink file prefix",type=str,required=True)
req_grp.add_argument("--tvec","-t",dest="stvec",help="input test vector file prefix",type=str,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
args=parser.parse_args()
print(args)

# Read in plink files
(bim_test, fam_test, G_test) = read_plink(args.inpre+"genos-test_common")
(bim_gwas, fam_gwas, G_gwas) = read_plink(args.inpre+"genos-gwas_common")

# Mean Center Test
X = G_test.transpose()
scaler = dask_ml.preprocessing.StandardScaler(with_mean=True, with_std=False).fit(X)
X = scaler.transform(X)

# Mean Center GWAS
M = G_gwas.transpose()
scaler = dask_ml.preprocessing.StandardScaler(with_mean=True, with_std=False).fit(M)
M = scaler.transform(M)

# Calculate Test covariance matrix 
test_cov = dask.array.cov(X)

# SVD of covariance matrix 
u, s, v = dask.array.linalg.svd(test_cov)
n = s.shape[0]

# Read in Test Vec
stvec = np.loadtxt(args.stvec)

# Make pseudoinverse 
vals = dask.array.diag(1/s[0:(n-1)])
pinv = np.dot(u[0:n,0:(n-1)], vals).dot(u[0:n,0:(n-1)].transpose()).compute()

# Calculate Tm
K = np.dot(M, X.transpose()) / (M.shape[1] - 1)
Tm = np.dot(K,pinv).dot(stvec)
Tm = Tm.compute()

# Write to file 
np.savetxt(args.outpre,Tm)

