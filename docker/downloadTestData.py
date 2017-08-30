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
	parser.add_argument('sc',help="subchallenge",choices=['sc1','sc2','sc3'])
	parser.add_argument("--express", help="Express lane", action="store_true", default=False)
	args = parser.parse_args()
	syn = synapseclient.login()
	downloadDir = os.path.dirname(os.path.abspath(__file__))
	testDataDir = os.path.join(downloadDir,"evaluation_data")
	
	if not os.path.exists(testDataDir):
		os.mkdir(testDataDir)
	os.system("rm -f %s/*" % testDataDir)
	if args.sc == 'sc1':
		downloadData(syn, "syn10164401",testDataDir)
	elif args.sc == 'sc2':
		downloadData(syn, "syn10139559",testDataDir)
		#downloadData(syn, "syn10139560",testDataDir)
	else:
		downloadData(syn, "syn10139567",testDataDir)
		#downloadData(syn, "syn10139568",testDataDir)

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
