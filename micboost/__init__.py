# coding: utf-8



import numpy as np
import pandas as pd
import random as rd
import xgboost as xgb
from sklearn.metrics import mutual_info_score
import math as mth
import scipy.linalg as la
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics.cluster import normalized_mutual_info_score




# This function calculates the mutual information between two sets of data, namely X and y
# by default X is a binary variables, y is also a binary variable if ybool==True, else it is a continuous variable
def mutualinfoAuto(X, y, ybool, normal):
    if normal == True:
        return normalized_mutual_info_score(X, y)
    else:
        if ybool == False:
            try:
                return mutual_info_regression(X.values.reshape(-1, 1), y, 'auto')[0]
            except:
                return mutual_info_regression(X.reshape(-1, 1), y, 'auto')[0]
        else:
            return mutual_info_score(X, y, contingency=None)




# This is a feature selection method based on this paper (http://ieeexplore.ieee.org/document/4749258/?part=1)
# It selects the feature which has the high mutual_info with y but low mutual_info with already selected features
def GscoreSelection(unselectedfeature, selectedfeature, dataS):
    temp = -1
    besti = 0
    if len(selectedfeature) == 0:
        for i in unselectedfeature:
            Gscore = mutualinfoAuto(dataS.iloc[:, i], dataS.iloc[:, 0], False, False)
        if Gscore > temp:
            temp = Gscore
            besti = i
    else:
        for i in unselectedfeature:
            NIsum = 0.0
            for j in selectedfeature:
                NIsum += mutualinfoAuto(dataS.iloc[:, i], dataS.iloc[:, j], True, True)
            Gscore = mutualinfoAuto(dataS.iloc[:, i], dataS.iloc[:, 0], False, False) - NIsum / len(selectedfeature)
            if Gscore > temp:
                temp = Gscore
                besti = i
    print besti
    print Gscore
    return besti


# besti is the index of the feature, it is an integer not a string




# This function uses GscoreSelecton to iterate over all features and select them based on Gscores
def featureSelection(f, dataS):
    s = 1

    unselectedfeature = range(1, len(dataS.columns) - 1, 1)
    selectedfeature = []
    while s < f:
        s += 1
        newfeature = GscoreSelection(unselectedfeature, selectedfeature, dataS)
        unselectedfeature.remove(newfeature)
        selectedfeature += [newfeature]

    return selectedfeature




# This function simulates the a table including strains with random snps and associated snps
# strainnum is the intended strain number, snp1, snp2 are the variance of snp1 and snp2, snp12 is the covariance between snp1 and snp2
# By default, the covariance between snp1 and snp2 is the average of their variance
def SNPSimulator(strainnum, snp1, snp2, snp12, randsnps):
    if snp12 is None:
        snp12 = (snp1 + snp2) / 2.0
    headlist = np.random.normal(0, 1, 3 * strainnum)
    newheadlist = np.array(headlist).reshape(3, strainnum)
    covmtrx = np.array([1, snp1, snp2, snp1, 1, snp12, snp2, snp12, 1]).reshape(3, 3)
    try:
        V = la.cholesky(covmtrx)
        cdgSNP12 = np.dot(V, newheadlist)
        for i in range(0, len(cdgSNP12[0]), 1):
            for j in range(1, 3, 1):
                if cdgSNP12[j][i] > 0:
                    cdgSNP12[j][i] = 1

                else:
                    cdgSNP12[j][i] = 0
        randomsnps = np.random.binomial(1, 0.5, size=strainnum * randsnps).reshape(randsnps, strainnum)
        result = np.row_stack((cdgSNP12, randomsnps))
        return np.transpose(result)
    except:
        return 'Please enter a valid covariance'


# sample=pd.DataFrame(SNPSimulator(30,0.5,0.5,0.2,100))
# sample.head()




# This tuning function uses the SNPSimulator function with parameter (30,0.5,0.5,0.2,100)
# and it returns the sorted result of xgboost
def tuning(max_depth, eta, numround, gm):
    sample = SNPSimulator(30, 0.5, 0.5, 0.2, 100)
    X = np.transpose(sample[1:-1])
    y = sample[0]
    dtrain = xgb.DMatrix(X, label=y)
    param = {'max_depth': max_depth, 'eta': eta, 'silent': 1, 'gamma': gm}
    num_round = numround
    bst = xgb.train(param, dtrain, num_round)
    score = bst.get_score(importance_type='gain')
    result = sorted(score.iteritems(), key=lambda (k, v): (v, k), reverse=True)
    return result


# tuning(1, 0.01, 10, 0.1)




# These two functions are used to tune the parameters of the xgboost
# We hope xgboost is able to detect the correlated snps in our simulated data
# Everytime xgboost get feature0 and feature1, we reward it by giving 1 or 0.5 scores

def xgscore(a, b, c, d):
    score = 0.0
    xgresult = tuning(a, b, c, d)
    try:
        if ((xgresult[0][0] == 'f0') and (xgresult[1][0] == 'f1')) or (
            (xgresult[0][0] == 'f1') and (xgresult[1][0] == 'f0')):
            score += 1
    except:
        pass

    try:
        if (xgresult[0][0] == 'f1') or (xgresult[0][0] == 'f0'):
            score += 0.5
    except:
        pass

    try:
        if (xgresult[1][0] == 'f0') or (xgresult[1][0] == 'f1'):
            score += 0.5
    except:
        pass
    return score


# Whichever combination of the parameters obtains the highest scores is considered the best combination,
# so we return this set of parameters
# This algo is kind of complicated and can be certainly improved
def getBestParameters():
    maxscore = 0
    parameterlist = []
    for depth in range(1, 10, 1):
        for eta in range(0, 10, 1):
            for numberround in range(0, 100, 5):
                for gm in range(0, 5, 1):
                    temp = 0.0
                    for time in range(0, 10, 1):
                        temp += xgscore(depth, float(eta) / 10.0, numberround, gm)
                        if temp > maxscore:
                            maxscore = temp
                            parameterlist = [depth, float(eta) / 10.0, numberround, gm]
    return parameterlist, maxscore


# getBestParameters()


#

# This function calculates the mutual information of a binary dataset with NA values
# The way it deals with NA values is to ignore the row with NA's in both snps
# and calculate the mutualinfo based on the rest of the value
def mutualinfo_dataframe(mydataframe):
    r, c = mydataframe.shape
    df1 = pd.DataFrame(mydataframe).copy()
    df1.iloc[:, :] = 0
    df = df1.corr()
    l = list(df.columns)
    tet25 = pd.read_csv('x-tet25.csv')
    for snp1 in l:
        ix = 0
        for snp2 in l:
            X = []
            y = []
            for i in range(0, r, 1):
                if (mydataframe[snp1][i] == mydataframe[snp1][i]) & (mydataframe[snp2][i] == mydataframe[snp2][i]):
                    X += [int(mydataframe[snp1][i])]
                    y += [int(mydataframe[snp2][i])]
            df[snp1][ix] = mutualinfoAuto(X, y, True, True)
            ix += 1
    return df


# tet25=pd.read_csv('x-tet25.csv')
# mutualinfo_dataframe(tet25.iloc[:,2:])




# This function is used to reform the mutualinfo matrix so that we can see which two snps have the highest mutual info
# or lowest mutual info
def reForm(mydataframe):
    muinfo = mydataframe
    r, c = mydataframe.shape
    names = list(mydataframe.columns)
    l = []
    for n in names:
        l += [n]
    muinfo["Unnamed: 0"] = l
    ar = np.array(muinfo["Unnamed: 0"])

    col = []
    for i in ar:
        col += [str(i)]
    col
    infoTable = pd.melt(muinfo, id_vars=['Unnamed: 0'], value_vars=col)
    infoTable = infoTable.rename(columns={'Unnamed: 0': 'SNPA', 'variable': 'SNPB'})
    infoTable['X'] = 0
    infoTable['Y'] = 0
    listX = []
    listY = []
    for i in range(0, r, 1):
        for j in range(0, r, 1):
            listX += [i + 1]
            listY += [j + 1]
    infoTable.loc[:, 'X'] = pd.Series(listX, index=infoTable.index)
    infoTable.loc[:, 'Y'] = pd.Series(listY, index=infoTable.index)
    return infoTable


# tet25=pd.read_csv('x-tet25.csv')
# tet25mutualInfo=mutualinfo_dataframe(tet25.iloc[:,2:])
# infoTable=reForm(tet25mutualInfo)
# uniqueTable=infoTable[infoTable["X"]!=infoTable["Y"]]
# uniqueTable=uniqueTable.sort_values(by='value', ascending=False)
# uniqueTable




# This function takes all the parameters of xgboost, and a input dataframe including the name of X column and y column.
# It uses cross validation to evaluate the result from xgboost and returns the test score as well as the significant snps.
def xgb_run(max_depth, eta, numround, gm, inputsample, inputfeature, colname):
    sample = inputsample
    X = sample[inputfeature]
    y = sample[colname]
    dtrain = xgb.DMatrix(X, label=y)
    param = {'max_depth': max_depth, 'eta': eta, 'silent': 1, 'gamma': gm}
    num_round = numround
    score = xgb.train(param, dtrain, num_round).get_score(importance_type='gain')
    return sorted(score.iteritems(), key=lambda (k, v): (v, k), reverse=True)


# df=pd.read_csv('x-azt32.csv')
# xgb_run(3,0.03, 3, .3,df,['SNP2','SNP3','SNP4'],'Azt32')




# This function takes all the parameters of xgboost, and a input dataframe including the name of X column and y column.
# It uses cross validation to evaluate the result from xgboost and returns the test score as well as the significant snps.
def xgtest(max_depth, eta, numround, gm, inputsample, inputfeature, colname):
    sample = inputsample
    X = sample[inputfeature]
    y = sample[colname]
    dtrain = xgb.DMatrix(X, label=y)
    param = {'max_depth': max_depth, 'eta': eta, 'silent': 1, 'gamma': gm}
    num_round = numround
    bst = xgb.cv(param, dtrain, num_round)
    score = xgb.train(param, dtrain, num_round).get_score(importance_type='gain')
    return bst



# This function takes a csvfile (with first column strainID, second column antibiotic resistence level,
# third column and after snp states) and the name of the second column and the mutualinfo_table of snps,
# and a threshold(rt) which is used to exclude high mutual info snps
# The porpose of this function is to reduce the dimension of the dataset by excluding all the high mutual info snps
# It initializes by select a snps which is mostly related to the y value and keep adding snps which has low mutual info
# with all already selected snps
# From the result we can see, the reduction in dimension slightly increases the xgboost test scores (given correct parameters)
def reductiontest(rt, csvname, colname, mutualinfo_table):
    df = pd.read_csv(csvname)
    names = df.columns.values[2:]
    l = []
    r, c = df.shape
    for n in names:
        l += [n]

    unselectedfeature = l
    selectedfeature = []
    temp = -1
    for j in unselectedfeature:
        X = []
        y = []
        for i in range(0, r, 1):
            if (df[j][i] == df[j][i]):
                X += [float(df.iloc[i][1])]
                y += [float(df[j][i])]
        X = np.array(X)
        score = mutualinfoAuto(X, y, False, False)
        if score > temp:
            temp = score
            besti = j
    firstfeature = j
    selectedfeature += [firstfeature]
    unselectedfeature.remove(firstfeature)
    for i in range(0, len(unselectedfeature), 1):
        newfeature = unselectedfeature[i]
        addflag = 1
        for feature in selectedfeature:
            if mutualinfo_table[feature][i + 1] > rt:
                addflag = 0
        if addflag == 1:
            selectedfeature += [newfeature]
    df = pd.read_csv(csvname)
    cdglevel = [colname]
    wholeselected = cdglevel + selectedfeature
    selected_df = df[wholeselected]
    resultdf = xgtest(3, 0.03, 3, .3, selected_df, selectedfeature, colname)
    return resultdf

# for i in [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]:
#     print i,reductiontest(i,'x-azt32.csv','Azt32',tet25mutualInfo)

