import csv
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw, MACCSkeys, rdFingerprintGenerator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

###User inputs###
csv_file = input('Provide a file path to a .csv file with two columns, ID and SMILES: ').strip().strip('"')
fingerprints = input ('What fingerprinting method should be used? MACCS or Morgan ').lower()
similarity_choice = input ('What similarity coefficient should be used? Tanimoto or Dice similarity ').lower()
distance_cutoff = float(input ('What distance cut off should be used? This is a measure of similarity, with 0 being identical, and 1 being no similarity '))
distribution = input('generate a bar chart showing number of molecules in each cluster? y/n ').lower()
visualise = input('Create a figure of the chemical structures in the largest cluster(s)? y/n ').lower()
excel = input('Create an .xlsx spead sheet of clusters? y/n ').lower()


###Extracting compounds###
with open(csv_file, newline="") as file:
	structures_file = csv.DictReader(file)

	fp_list = []
	
	for col in structures_file.fieldnames:
		if 'ID' in col.upper():
			id_column = col
			break

	for col in structures_file.fieldnames:
		if 'SMILES' in col.upper():
			SMILES_column = col
			break

	for row in structures_file:
		title = row[id_column]
		
		smiles = row[SMILES_column]

		###Generating fingerprints, MACCS, Morgan, or repeating if input was wrong###
		mol_gen = Chem.MolFromSmiles(smiles)

		if mol_gen is not None:

			while True:
			
				if fingerprints == 'maccs':
					maccs_fp = MACCSkeys.GenMACCSKeys(mol_gen)
					fp_list.append([title, smiles, maccs_fp])
					break

				elif fingerprints =='morgan':
					fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
					morgan_fp = fpg.GetCountFingerprint(mol_gen)
					fp_list.append([title, smiles, morgan_fp])
					break

				else:
					fingerprints = input ('Fingerprinting method not recognised, what fingerprinting method should be used? MACCS or Morgan ').lower()

	###Calculating similarity between every compound###
	distance_matrix = []
	
	for i in range(0, (len(fp_list))):
		for j in range(i):
			fp_1 = fp_list[i][2]
			fp_2 = fp_list[j][2]

			###Selecting similarity coefficient###
			while True:

				if similarity_choice == 'tanimoto':
					similarity = DataStructs.TanimotoSimilarity(fp_1, fp_2)
					break

				elif similarity_choice == 'dice':
					similarity = DataStructs.DiceSimilarity(fp_1, fp_2)
					break

				else:
					similarity_choice = input ('similarity coefficient not recognised, what similarity coefficient should be used? Tanimoto or Dice similarity ').lower()


			distance = 1 - similarity
			distance_matrix.append(distance)

	###Making sure the distance cutoff is between 0 and 1###

	while True:
		if 0 <= distance_cutoff <=1:
			break
		else:
			distance_cutoff = float(input ('distance cut off outside allowed range, what distance cut off should be used? This is a measure of similarity, with 0 being identical, and 1 being no similarity '))


	###Butina clustering###
	clusters = Butina.ClusterData(distance_matrix, len(fp_list), distance_cutoff, isDistData=True)
	no_clusters = len(clusters)
	clusters = sorted(clusters, key=len, reverse=True)

	###Creating distribution plot###
	if distribution == 'y':
		fig, ax = plt.subplots(figsize=(15, 4))
		ax.set_title(fingerprints + " " + similarity_choice + " " + str(distance_cutoff))
		ax.set_xlabel("Cluster index")
		ax.set_ylabel("Number of molecules")
		ax.set_xlim(0, len(clusters))
		ax.set_xticks(range(len(clusters))) 
		ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
		ax.set_yticks(range(len(clusters[0])))
		ax.yaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
		ax.bar(range(0, len(clusters)), [len(cluster) for cluster in clusters], lw=5)
		plt.show()

	all_dicts = {}

	for i in range(0, no_clusters):

		for compound in clusters[i]:
			molecule = fp_list[compound][0]
			molecule_smiles = fp_list[compound][1]
			molecule_mol_gen = Chem.MolFromSmiles(molecule_smiles)
			if i not in all_dicts:
				all_dicts[i]={molecule: molecule_mol_gen}

			else:
				all_dicts[i][molecule] = molecule_mol_gen

	rows=[]

	for key, value in all_dicts.items():
		row=[key]
		for k, v in value.items():
			row.extend([k, v])

		rows.append(row)

	###creating figures###

	if visualise == 'y':

		print('There are ' + str(no_clusters) + ' clusters')
		no_clusters_visualise = int(float(input('How many clusters should be visualised in the figure? Note clusters will be visualised from largest to smallest ')))

		while True:
			if 0 < no_clusters_visualise <= no_clusters:
				break

			else:
				no_clusters_visualise = int(float(input('How many clusters should be visualised? Must be smaller than the total number of clusters ' + str(no_clusters) + ' ')))

		for i in range(0, no_clusters_visualise):
 			img=Draw.MolsToGridImage([rows[i][j] for j in range(2, len(rows[i]), 2)], molsPerRow=10,  subImgSize=(200,200), legends=[rows[i][j] for j in range(1, len(rows[i]), 2)])
 			img.show()

 	###creating an .xlsx file###

	if excel == 'y':

		name = fingerprints + "_" + similarity_choice + ".xlsx"

		excel = []

		for cluster in rows:
			excel.append([cluster[i] for i in range(1, len(cluster), 2)])

		df = pd.DataFrame(excel)
		df.to_excel(name, index=False, header=False)