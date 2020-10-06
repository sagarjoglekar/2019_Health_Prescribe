import pandas as pd 
import glob
import csv
import os
import tarfile
import inspect
import re
import numpy as np


#Different functions to extract different data form the prescription 

def cleanStringofUTF(string):
    cleaned = string.replace('\xe8','e').replace('\xf6','o')
    return cleaned

def enrichdrugs(chem_dict , drugs):
    diabetes_drug_words = [drugs[k]['name'].lower() for k in drugs]
    for drug in chem_dict:
        Name = chem_dict[drug]['name'].replace('(','').replace(')','')
        slot1 = Name.lower().split('/')
        slot2 = Name.lower().split(' ')
        slot3 = Name.lower().split(' & ')
        common1 = set(diabetes_drug_words).intersection(slot1)
        common2 = set(diabetes_drug_words).intersection(slot2)
        common3 = set(diabetes_drug_words).intersection(slot3)
        
        if len(common1) > 0 or len(common2) > 0 or len(common3) > 0:
            print(common1 , common2 , common3)
            drugs[chem_dict[drug]['code']] = {'disease':'' , 'disease_given_drug':0.0 , 'matched_disease':'', 'name':chem_dict[drug]['name'].strip() }

            
            
def makeChemDict(BNF_Chem):
    chem_dict = {}
    for index, row in BNF_Chem.iterrows():
        chem_dict[row['UNII_drugbank']] = {}
        chem_dict[row['UNII_drugbank']]['name'] = row['NAME']
        chem_dict[row['UNII_drugbank']]['code'] = row['BNF_code']
    return chem_dict
    
def getDrugCategory(categorylist, BNF_Chem, drugbankDict):
    allMatched = []
    drugs = {}
    chem_dict = makeChemDict(BNF_Chem)
    
    for k in drugbankDict:
        if len(drugbankDict[k]['Categories']) > 0:
            for cat in drugbankDict[k]['Categories']:
                matched_memo = []
                catString = cat.values()[0]#.split('\u2014')[-1]
                t = catString.lower().strip()
                for categoryString in categorylist:
                    categoryString = categoryString.lower()
                    if t.find(categoryString) >= 0:
                        matched_memo.append(categoryString)
                if k in chem_dict:
                    if len(matched_memo) > 0:# == len(categorylist):
                        allMatched.append(k)
                        print(chem_dict[k])
                        drugs[chem_dict[k]['code']] = {}
                        drugs[chem_dict[k]['code']]['name'] = chem_dict[k]['name']
                        drugs[chem_dict[k]['code']]['matched_cat'] = categorylist
    enrichdrugs(chem_dict,drugs)               
    return list(set(allMatched)) , drugs


def getDrugforDiseaseDrugbank(categorylist, BNF_Chem, drugbankDict):
    allMatched = []
    drugs = {}
    chem_dict = makeChemDict(BNF_Chem)
    
    for k in drugbankDict:
        if len(drugbankDict[k]['Associations']) > 0:
            for cat in drugbankDict[k]['Associations']:
                matched_memo = []
                catString = cat.values()[0]
                t = catString.lower().strip()
                for categoryString in categorylist:
                    categoryString = categoryString.lower()
                    if t.find(categoryString) >= 0:
                        matched_memo.append(categoryString)
                if k in chem_dict:
                    if len(matched_memo) > 0:
                        allMatched.append(k)
                        print(chem_dict[k])
                        drugs[chem_dict[k]['code']] = {}
                        drugs[chem_dict[k]['code']]['name'] = chem_dict[k]['name']
                        drugs[chem_dict[k]['code']]['matched_cat'] = categorylist
    enrichdrugs(chem_dict,drugs)               
    return  allMatched , drugs


def findDrugsForDisease(Graph, Disease, BNF_Chem ):#,threshProb):
    chem_dict = makeChemDict(BNF_Chem)
    drugs = {}
    for e in Graph.edges(data=True):
        if (cleanStringofUTF(e[1]).lower().find(Disease.lower()) >=0) or (cleanStringofUTF(e[0]).lower().find(Disease.lower()) >= 0) :
            drugNode = ''
            matchedDisease = ''
            if Graph.node[e[0]]['type'] == 'symptom':
                drugNode = e[1]
                matchedDisease = e[0]
            else:
                drugNode = e[0]
                matchedDisease = e[1]
            drugs[Graph.node[drugNode]['Id']] = {}
            drugs[Graph.node[drugNode]['Id']]['name'] = drugNode
            drugs[Graph.node[drugNode]['Id']]['matched_disease'] = matchedDisease
            drugs[Graph.node[drugNode]['Id']]['disease'] = Disease
    enrichdrugs(chem_dict,drugs)
    return drugs




def findDrugsForCategory(Graph, Cat, BNF_Chem ):#,threshProb):
    chem_dict = makeChemDict(BNF_Chem)
    drugs = {}
    for e in Graph.edges(data=True):
        if (cleanStringofUTF(e[1]).lower().find(Cat.lower()) >=0) or (cleanStringofUTF(e[0]).lower().find(Cat.lower()) >= 0) :
            drugNode = ''
            matchedDisease = ''
            if Graph.node[e[0]]['type'] == 'category':
                drugNode = e[1]
                matchedDisease = e[0]
            else:
                drugNode = e[0]
                matchedDisease = e[1]
            print(Graph.node[drugNode]['Id'])
            drugs[Graph.node[drugNode]['Id']] = {}
            drugs[Graph.node[drugNode]['Id']]['name'] = drugNode
            drugs[Graph.node[drugNode]['Id']]['matched_cat'] = matchedDisease
            drugs[Graph.node[drugNode]['Id']]['category'] = Cat
    enrichdrugs(chem_dict,drugs)
    return drugs

def findDrugByName(Graph, name):
    foundDrugs = []
    for node in Graph.nodes(data=True):
        if node[1]['type'] == 'drug':
            if name.lower() in node[0].lower().split(' '):
                foundDrugs.append(node)
    return foundDrugs

def relabelNodes(graph):
    mapping = {}
    for node in graph.nodes():
        mapping[node] = cleanStringofUTF(node)
        
    return nx.relabel_nodes(graph,mapping)


def readChem(year,files):
    for f in files:
        if f.find(year)>0:
            if f.find('CHEM') > 0:
                df = pd.read_csv(f,compression='infer',header=None,skiprows=1)
                return df

def readAddress(year,files):
    for f in files:
        if f.find(year)>0:
            if f.find('ADDR') > 0:
                df = pd.read_csv(f,compression='infer',header=None,skiprows=1)
                return df

def readPrescriptions(year,files):
    for f in files:
        if f.find(year)>0:
            if f.find('PDPI') > 0:
                df = pd.read_csv(f,compression='infer',header=None,skiprows=1)
                return df

def extractPostCodesDict(addrDf):
    postcodeDict = {}
    for index,row in addrDf.iterrows():
        try:
            postcodeDict[row[1]] = row[7].strip()
        except:
            print("Found invalid entry")
    return postcodeDict

def checkIndex(index):
    if index%100 == 0:
        return True
    else:
        return False

def getPC(key , postcodeDict):
    codes = []
    for k in key:
        if k in postcodeDict:
            codes.append(postcodeDict[k])
        else:
            codes.append('')
    return pd.Series(codes,index=key.index)

def getPostcode(df,postcodeDict):
    df[10] = ''
    df[10] = df.groupby(2)[2].apply(getPC,postcodeDict)   
    return df

def getDrugFamily(key, diseaseMap):
    found = 'N/A'
    for dcode in diseaseMap: 
        if key.name.find(dcode) == 0:
            found = dcode
            break
    drug = [found]*len(key)
    return pd.Series(drug,index=key.index)

def getDisease(key, diseaseMap):
    found = 'N/A'
    for dcode in diseaseMap: 
        if key.name.find(dcode) == 0:
            found = diseaseMap[dcode]['disease'].replace('\"','').replace('+',' ')
            break
    drug = [found]*len(key)
    return pd.Series(drug,index=key.index)

def getDrug(key, diseaseMap):
    found = 'N/A'
    for dcode in diseaseMap: 
        if key.name.find(dcode) == 0:
            found = diseaseMap[dcode]['chemName']
            break
    drug = [found]*len(key)
    return pd.Series(drug,index=key.index)

def getDrugPotency(key):
    name = list(set(key))
    if len(name) > 1:
        print("found synonyms")
    text= name[-1]
    found = 0.0
    switch1 = text.find('mcg')
    switch2 = text.find('mg')
    switch3 = text.find('ml')
    
    if switch1 > 0 or switch2 > 0 or switch3 > 0:
        weight = re.findall(r'[0-9]*\.?[0-9]+', text)
        if len(weight) > 0:
            found = max([float(k) for k in weight])
            if switch1 > 0:
                found = found/1000.0
    potency = [found]*len(key)
    return pd.Series(potency,index=key.index)

def countSpecificDrugs(Df, drugs,GPs):
    df_slice = Df.groupby(3)[3].apply(getDrug,drugs)
    selected = df_slice[df_slice!='N/A']
    len(selected)
    df_selected =  Df.iloc[selected.index,:]
    df_selected = df_selected[df_selected[2].isin(GPs.keys())]
    return np.sum(df_selected[5])


def countSpecificDrugCosts(Df, drugs,GPs):
    df_slice = Df.groupby(3)[3].apply(getDrug,drugs)
    selected = df_slice[df_slice!='N/A']
    len(selected)
    df_selected =  Df.iloc[selected.index,:]
    df_selected = df_selected[df_selected[2].isin(GPs.keys())]
    return np.sum(df_selected[6])

def countDrugsByCategoryList(pdp,codes):
    total_drugs = 0.0
    drugs = None
    for name , group in pdp.groupby('3'):
        for dcode in codes:
            if name.find(dcode) == 0:
                total_drugs+=np.sum(group['5'])
    return total_drugs

def countPrescriptionsByCategoryList(pdp,codes):
    total_drugs = 0.0
    drugs = None
    for name , group in pdp.groupby('3'):
        for dcode in codes:
            if name.find(dcode) == 0:
                total_drugs+=len(group)
    return total_drugs


def countDrugsCostByCategoryList(pdp,codes):
    total_drugs = 0.0
    for name , group in pdp.groupby('3'):
        for dcode in codes:
            if name.find(dcode) == 0:
                total_drugs+=np.sum(group['7'])
    return total_drugs

def countDrugsCostByGenerics(pdp,codes):
    total_drugs = 0.0
    generics= pdp[pdp['20'] == 'AA']
    for name , group in generics.groupby('3'):
        for dcode in codes:
            if name.find(dcode) == 0:
                total_drugs+=np.sum(group['7'])
    return total_drugs

def compareCostsForGenericsAndBranded(pdp,codes):
    genericsCosts = {}
    brandedCosts = {}
    for name , group in pdp.groupby('16'):
        for dcode in codes:
            if name == dcode:
                generics= group[group['20'] == 'AA']
                if len(generics)>0:
                    total_drugs =np.sum(group['7'])
                    generic_drugs = np.sum(generics['7'])
                    brandedCosts[dcode] = total_drugs - generic_drugs
                    genericsCosts[dcode] = generic_drugs
                else:
                    print("Did not find any generic drugs")
                    continue
    return brandedCosts , genericsCosts


def countTotalDrugDosage(pdp,codes):
    total_drugs = 0.0
    for name , group in pdp.groupby('3'):
        for dcode in codes:
            if name.find(dcode) == 0:
                total_drugs+=np.sum(group['19'])
    return total_drugs


def loadGraphs():
    
    with open('Bipartite_Drug_graph.pkl','rb') as f:
        drug_association_graph  = pkl.load(f)

    with open('Bipartite_Drug_category_graph.pkl','rb') as f:
        drug_cat_association_graph  = pkl.load(f)
        


            
