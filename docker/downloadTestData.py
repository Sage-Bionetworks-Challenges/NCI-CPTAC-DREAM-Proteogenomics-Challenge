import synapseclient
import argparse
import os
import shutil

def main():
	parser = argparse.ArgumentParser()
	parser.add_arguemnt('sc',help="subchallenge",choices=['sc1','sc2','sc3'])
	parser.add_argument("--express", help="Express lane", action="store_true", default=False)
	args = parser.parse_args()
	syn = synapseclient.login()
	testDataDir = os.path.join(os.getenv("HOME"),"evaluation_data")
	
	if not os.path.exists(testDataDir):
		os.mkdir(testDataDir)
	# cna = syn.get("syn10067991")
	# if args.sc in ['sc1', 'sc2']:
	# 	proteome = syn.get("syn10134684")
	# else:
	# 	proteome = syn.get("syn10067993")
	# rna = syn.get("syn10067992")

	# if args.express:
	# 	shutil.copy(cna.path, testDataDir)
	# 	shutil.copy(proteome.path, "%s/pros_ova_proteome_sort_all_gene.txt" % testDataDir)
	# 	shutil.copy(rna.path, testDataDir)
	# else:
	# 	shutil.copy(cna.path, testDataDir)
	# 	shutil.copy(proteome.path,  "%s/pros_ova_proteome_sort_all_gene.txt" % testDataDir)
	# 	shutil.copy(rna.path, testDataDir)

	cna = syn.get("syn9974010")
	if args.sc in ['sc1', 'sc2']:
		proteome = syn.get("syn10134685")
	else:
		proteome = syn.get("syn9974009")
	rna = syn.get("syn9974005")

	if args.express:
		shutil.copy(cna.path, testDataDir)
		shutil.copy(proteome.path, "%s/pros_ova_proteome_sort_common_gene_6577.txt" % testDataDir)
		shutil.copy(rna.path, testDataDir)
	else:
		shutil.copy(cna.path, testDataDir)
		shutil.copy(proteome.path,  "%s/pros_ova_proteome_sort_common_gene_6577.txt" % testDataDir)
		shutil.copy(rna.path, testDataDir)


if __name__ == '__main__':
	main()
