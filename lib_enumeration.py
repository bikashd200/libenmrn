#!/usr/bin/env python3
import os
import sys
from rdkit import Chem
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import pandas as pd
import csv
#from rdkit import Chem
from rdkit.Chem import rdChemReactions



params = Chem.SmilesParserParams()
params.removeHs=False # draw and work with explicit Hs

suppl1 = Chem.SmilesMolSupplier(sys.argv[1])#'covalent_compounds_with_GSH_reactivity.smi'
suppl2 = Chem.SmilesMolSupplier(sys.argv[2])

def matchCompound(subele, supplele):
    matches = supplele.GetSubstructMatches(Chem.RemoveHs(subele))
    #print(matches)
    if matches == ():
        return
    atomset = list(sum(list(matches), ()))
    m1 = Chem.AddHs(supplele, onlyOnAtoms=atomset)
    AllChem.Compute2DCoords(m1)
    #print(m1)
    return m1



def atomMapping(m_1, m1, ms, allsubs):
    sub = m1.GetSubstructMatches(m_1)
    if len(sub) > 0:
        ms.append(m1)
        allsubs.append(sub[0])

#m = Chem.MolFromMolFile('Cys_3D.mol', removeHs=False)
prod = Chem.MolFromMolFile('thioether_prod.mol', removeHs=True)
#m = Chem.MolFromMolFile('glutathione.mol', removeHs=False)

rxn1 = rdChemReactions.ReactionFromSmarts('[S:1][H].[N:2][C:3](=[O:4])[C:5](=[C:6])>>[N:2][C:3](=[O:4])[C:5]([H])[C:6][S:1]') #NC(=O)C=C



import csv
file_count=0
#count=0
subs_map = []
final_ps =[]
#usmis = set()
usmis = set()
#with open(sys.argv[2], 'w') as f:
    #'Cys_acrylamides.smi'
#    out_file = csv.writer(f)
for m1 in suppl1:
    if m1 is None:
        continue
    #usmis = set()
    ID1 = m1.GetProp("_Name")
    print(ID1)
    m1 = Chem.AddHs(m1)
    for m2 in suppl2:
        if m2 is None:
            continue
        ID2 = m2.GetProp("_Name")
        print(ID2)
        m2 = Chem.AddHs(m2)
        products = rxn1.RunReactants((m1,m2))
        #count = 1
        for p in products:
            smi = Chem.MolToSmiles(Chem.RemoveHs(p[0]))
            print(smi)
            #smile_string = ([smi + " " + ID1 + '_' + ID2 ])
            smile_string = smi + " " + ID1 + '_' + ID2
            #count+=1
            mol_p = p[0]
            if mol_p is None:
                continue
            usmis.add(smile_string)

for i in usmis:
    mol_p = Chem.MolFromSmiles(i)
    #display(mol_p)
    #count+=1
    #mol_p.SetProp("ID",'Prod_'+ str(count))
    #writer.write(mol_p)
    atomMapping(prod, mol_p, final_ps, subs_map)

count = 0
print(len(final_ps))
with open(sys.argv[3], 'w') as f:
    out_file = csv.writer(f)
    for molp in final_ps:
        smi = Chem.MolToSmiles(molp)
        count+=1       
        smile_string = ([smi + ' ' + 'prod_'+ str(count)])
        #f.write(smile_string)
          
        out_file.writerow(smile_string)
        #count+=1       
f.close() 
