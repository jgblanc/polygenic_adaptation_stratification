# This script generates Test vector and calculates Tm 
print("I am running!")

# Import modules
import argparse
from pandas_plink import read_plink
import pandas as pd
import dask
import dask_ml.preprocessing
import numpy as np
import dask.array as da

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
print("Read in Plink Files")

# Mean Center Test
X = G_test.transpose()
scaler = dask_ml.preprocessing.StandardScaler(with_mean=True, with_std=False).fit(X)
X = scaler.transform(X)
print("Mean centered X")

# Mean Center GWAS
M = G_gwas.transpose()
scaler = dask_ml.preprocessing.StandardScaler(with_mean=True, with_std=False).fit(M)
M = scaler.transform(M)
print("Mean centered M") 

# Calculate Test covariance matrix 
test_cov = dask.array.cov(X)
print("MAde covariance matrix")

# SVD of covariance matrix 
u, s, v = dask.array.linalg.svd_compressed(test_cov, test_cov.shape[0])
n = s.shape[0]
print("Did SVD")

# Read in Test Vec
stvec = np.loadtxt(args.stvec)
print("Read Test Vec")

# Make pseudoinverse 
vals = dask.array.diag(1/s[0:(n-1)])
pinv = np.dot(u[0:n,0:(n-1)], vals).dot(u[0:n,0:(n-1)].transpose())
print("Made pseuodoinverse")

# Calculate Tm
K = np.dot(M, X.transpose()) / (M.shape[1] - 1)
Tm = np.dot(K,pinv).dot(stvec)
print("Calculated Tm")

# Write to file 
da.to_hdf5(args.outpre, '/output', Tm)
print("Wrote hdf5 file")

#Tm = Tm.compute()
#np.savetxt("test.txt",Tm)
#file_name=args.outpre
#outF = open(file_name, "w")
#for i in range(Tm.shape[0]):
#    outF.write(str(Tm[i].compute()))
#outF.close()
