#!/usr/bin/python
'''USAGE: splitMetab.py
This script splits metabolomic data and creates python dictionaries for each condition
'''
import sys
import pickle
import numpy

def geo_mean(values):
    x = numpy.array(values)
    x = x.prod() ** (1.0/len(x))
    return x


biochemical_set = set()

strep_biochem_dict = {}
cef_biochem_dict = {}
clinda_biochem_dict = {}
gf_biochem_dict = {}

strep_pubchem_dict = {}
cef_pubchem_dict = {}
clinda_pubchem_dict = {}
gf_pubchem_dict = {}


strep_kegg_dict = {}
cef_kegg_dict = {}
clinda_kegg_dict = {}
gf_kegg_dict = {}


with open(sys.argv[1], 'r') as metabolomes:

	for entry in metabolomes:
		entry = entry.split()

		if entry[0] == 'BIOCHEMICAL':
			continue

		biochemical = str(entry[0])
		pubchem = str(entry[3])
		kegg = str(entry[4])

		cef_entry = geo_mean([float(x) for x in entry[14:23]])
		clinda_entry = geo_mean([float(x) for x in entry[32:41]])
		strep_entry = geo_mean([float(x) for x in entry[50:59]])
		gf_entry = geo_mean([float(x) for x in entry[68:77]])

		if not biochemical in biochemical_set:
			biochemical_set.add(biochemical)

			cef_biochem_dict[biochemical] = cef_entry
			clinda_biochem_dict[biochemical] = clinda_entry
			strep_biochem_dict[biochemical] = strep_entry
			gf_biochem_dict[biochemical] = gf_entry

			if pubchem != 'NA':
				cef_pubchem_dict[biochemical] = cef_entry
				clinda_pubchem_dict[biochemical] = clinda_entry
				strep_pubchem_dict[biochemical] = strep_entry
				gf_pubchem_dict[biochemical] = gf_entry
			if kegg != 'NA':
				cef_kegg_dict[biochemical] = cef_entry
				clinda_kegg_dict[biochemical] = clinda_entry
				strep_kegg_dict[biochemical] = strep_entry
				gf_kegg_dict[biochemical] = gf_entry

		else:
			cef_biochem_dict[biochemical] = cef_entry + cef_biochem_dict[biochemical]
			clinda_biochem_dict[biochemical] = clinda_entry + clinda_biochem_dict[biochemical]
			strep_biochem_dict[biochemical] = strep_entry + strep_biochem_dict[biochemical]
			gf_biochem_dict[biochemical] = gf_entry + gf_biochem_dict[biochemical]

			if pubchem != 'NA':
				cef_pubchem_dict[biochemical] = cef_entry + cef_pubchem_dict[biochemical]
				clinda_pubchem_dict[biochemical] = clinda_entry + clinda_pubchem_dict[biochemical]
				strep_pubchem_dict[biochemical] = strep_entry + strep_pubchem_dict[biochemical]
				gf_pubchem_dict[biochemical] = gf_entry + gf_pubchem_dict[biochemical]
			if kegg != 'NA':
				cef_kegg_dict[biochemical] = cef_entry + cef_kegg_dict[biochemical]
				clinda_kegg_dict[biochemical] = clinda_entry + clinda_kegg_dict[biochemical]
				strep_kegg_dict[biochemical] = strep_entry + strep_kegg_dict[biochemical]
				gf_kegg_dict[biochemical] = gf_entry + gf_kegg_dict[biochemical]


pickle.dump(cef_biochem_dict, open('cef_biochem_dict.pickle', 'wb'))
pickle.dump(clinda_biochem_dict, open('clinda_biochem_dict.pickle', 'wb'))
pickle.dump(strep_biochem_dict, open('strep_biochem_dict.pickle', 'wb'))
pickle.dump(gf_biochem_dict, open('gf_biochem_dict.pickle', 'wb'))

pickle.dump(cef_pubchem_dict, open('cef_pubchem_dict.pickle', 'wb'))
pickle.dump(clinda_pubchem_dict, open('clinda_pubchem_dict.pickle', 'wb'))
pickle.dump(strep_pubchem_dict, open('strep_pubchem_dict.pickle', 'wb'))
pickle.dump(gf_pubchem_dict, open('gf_pubchem_dict.pickle', 'wb'))

pickle.dump(cef_kegg_dict, open('cef_kegg_dict.pickle', 'wb'))
pickle.dump(clinda_kegg_dict, open('clinda_kegg_dict.pickle', 'wb'))
pickle.dump(strep_kegg_dict, open('strep_kegg_dict.pickle', 'wb'))
pickle.dump(gf_kegg_dict, open('gf_kegg_dict.pickle', 'wb'))


