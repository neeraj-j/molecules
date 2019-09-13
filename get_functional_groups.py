# find functional group using rdkit
# single thread performance is better than multithreaded performance !?!?!? 
# This code is derived from https://gist.github.com/CKannas/8954290


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import FunctionalGroups

from time import time
import csv

start = time()
# open target csv file
csvfile = open("fg_smils.csv", 'a')
header = ["pubchemId","functionalGroup", "SMILES"]
writer = csv.writer(csvfile)
writer.writerow(header)

smilesfile = "CID-SMILES"

fgs = FunctionalGroups.BuildFuncGroupHierarchy()


#for filename in filenames:
    
suppl = Chem.SmilesMolSupplier(smilesfile, delimiter=",",
                               smilesColumn=1, nameColumn=0, titleLine=False )

mols = [x for x in suppl if x is not None]
del suppl

print("processing %s with %d valid compounds" % (filename, (len(mols))))

#Draw.MolsToGridImage(mols[:20], molsPerRow=4, legends=[x.GetProp('_Name') for x in mols])

'''
# Dont want prints
print("Factional Groups")
print("Label\tSMARTS")
print("===============")
for x in fgs:
    print(x.label, "\t", x.smarts)
    for y in x.children:
        print(y.label, "\t", y.smarts)
        for z in y.children:
            print(z.label, "\t", z.smarts)
'''

zbbAllFGs = {}

def getFGs(fgs, res):
    if not fgs:
        return
    for x in fgs:
        patt = x.pattern
        tmp = [m for m in mols if m.HasSubstructMatch(patt)]
        # if there are functional groups then check its children also
        if len(tmp):
            res[x.label] = {'mols': tmp, 'pattern': patt}
            #print(x.label, Chem.MolToSmarts(patt))
            getFGs(x.children, res)
        # end if
    # end for

    return


getFGs(fgs, zbbAllFGs)
len(zbbAllFGs)
#mols[0].GetProp('_Name')
totalFGs = 0
for fgName in sorted(zbbAllFGs.keys()):
    totalFGs += len(zbbAllFGs[fgName]['mols'])
    #print("%s: Found %d" %(fgName, len(zbbAllFGs[fgName]['mols'])))
    for ml in zbbAllFGs[fgName]['mols']:
        id = ml.GetProp('_Name')
        smiles = Chem.MolToSmiles(ml)
        row = [id, fgName, smiles]
        writer.writerow(row)
    csvfile.flush()


del mols
del zbbAllFGs
print("Total number of functional groups: %d" % (totalFGs))

csvfile.close()
end = time() - start
print(" TOtal time taken: ", end)
