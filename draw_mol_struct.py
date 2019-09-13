
# draw molecular structure from fg_smiles.csv file 
# This code is derived from https://gist.github.com/CKannas/8954290

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import FunctionalGroups

from time import time
import pandas as pd
import matplotlib.pyplot as plt
import os
start = time()

baseImgDir = "./images/"

fgs = FunctionalGroups.BuildFuncGroupHierarchy()
# make dictionary of template
def makeTemplateDic(fgs, tempDic):
    for x in fgs:
        tempDic[x.label] = x.pattern
        #print(x.label, Chem.MolToSmarts(patt))
        makeTemplateDic(x.children, tempDic)
        # end if
    # end for
    return

templateDic = {}

makeTemplateDic(fgs, templateDic)

pd_iterator = pd.read_csv("./fg_smils.csv",  chunksize=20000)

## Get number of elements in each functional group. Only one time
fg_count = {}
for rowsDf in pd_iterator:
    rows = rowsDf.functionalGroup.value_counts()
    for row in rows.iteritems():
        if not row[0] in fg_count:
            fg_count[row[0]] = row[1]
        else:
            fg_count[row[0]] += row[1]
print(fg_count)
############# End ################

usedfgs = ['Amine.Cyclic','Alcohol.Aliphatic','Halogen.Aromatic','CarboxylicAcid.AlphaAmino']
imgcount = {'Amine.Cyclic':0,'Alcohol.Aliphatic':0,'Halogen.Aromatic':0,'CarboxylicAcid.AlphaAmino':0}
for rowsDf in pd_iterator:
    fgs = rowsDf.functionalGroup.unique()
    for ufg in fgs:
        if not ufg in usedfgs:
            continue
        if imgcount[ufg] >= 50000: # dont create more than 50K images per fg
            continue
        AllChem.Compute2DCoords(templateDic[ufg])
        rows = rowsDf[rowsDf['functionalGroup'] == ufg]
        for _, row in rows.iterrows():
            pubId = row["pubchemId"]
            smiles = row["SMILES"]
            mol = Chem.MolFromSmiles(smiles)
            AllChem.GenerateDepictionMatching2DStructure(mol,templateDic[ufg])
            dstdir = baseImgDir + ufg
            if not os.path.exists(dstdir):
                os.mkdir(dstdir)
            dstfile = dstdir +'/' + str(pubId) + ".png"
            Draw.MolToFile(mol,dstfile, size=(128, 128))
            imgcount[ufg] +=1
            #img = Draw.MolToImage(mol, size=(128, 128))
            #plt.imshow(img)

