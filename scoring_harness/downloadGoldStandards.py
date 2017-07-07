import synapseclient
import argparse
import os
import shutil

def downloadData(syn, parentId, testDataDir):
	download = syn.getChildren(parentId)
	for i in download:
		temp = syn.get(i['id'])
		shutil.copy(temp.path, "%s/%s" % (testDataDir, os.path.basename(temp.path)))

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--express", help="Express lane", action="store_true", default=False)
	args = parser.parse_args()
	syn = synapseclient.login()
	downloadDir = os.path.dirname(os.path.abspath(__file__))
	
	downloadData(syn, "syn10165897",downloadDir)
	downloadData(syn, "syn10165898",downloadDir)
	downloadData(syn, "syn10165900",downloadDir)

	# if args.express:
	# 	shutil.copy(cna.path, testDataDir)
	# 	shutil.copy(proteome.path, "%s/pros_ova_proteome_sort_common_gene_6577.txt" % testDataDir)
	# 	shutil.copy(rna.path, testDataDir)
	# else:
	# 	shutil.copy(cna.path, testDataDir)
	# 	shutil.copy(proteome.path,  "%s/pros_ova_proteome_sort_common_gene_6577.txt" % testDataDir)
	# 	shutil.copy(rna.path, testDataDir)


if __name__ == '__main__':
	main()
