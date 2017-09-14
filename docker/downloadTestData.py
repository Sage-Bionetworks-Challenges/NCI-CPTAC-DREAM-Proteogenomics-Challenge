import synapseclient
import argparse
import os
import shutil

def downloadData(syn, parentId, testDataDir, replaceName=None):
	download = syn.getChildren(parentId)
	for i in download:
		temp = syn.get(i['id'])
		if replaceName is not None:
			newName = os.path.basename(temp.path).replace(replaceName,'')
		else:
			newName = os.path.basename(temp.path)
		shutil.copy(temp.path, "%s/%s" % (testDataDir, newName))

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('sc',help="subchallenge",choices=['sc1','sc2','sc3'])
	parser.add_argument('round',help="round",choices=['1','2','final','test'])
	parser.add_argument("--express", help="Express lane", action="store_true", default=False)
	args = parser.parse_args()
	syn = synapseclient.login()
	downloadDir = os.path.dirname(os.path.abspath(__file__))
	testDataDir = os.path.join(downloadDir,"evaluation_data")
	
	if not os.path.exists(testDataDir):
		os.mkdir(testDataDir)
	os.system("rm -f %s/*" % testDataDir)

	if args.round == '1':
		sc2 = "syn10617839"
		sc3 = 'syn10617842'
		replace = "_round_1"
	elif args.round == '2':
		sc2 = 'syn10617840'
		sc3 = 'syn10617843'
		replace = "_round_2"
	elif args.round == 'test':
		sc2 = 'syn10139559'
		sc3 = 'syn10139567'
	else:
		sc2 = 'syn10617841'
		sc3 = 'syn10617844'
		replace = "_final_round"

	if args.sc == 'sc1':
		downloadData(syn, "syn10164401",testDataDir)
	elif args.sc == 'sc2':
		#downloadData(syn, sc2, testDataDir,replace)
		downloadData(syn, sc2, testDataDir)
	else:
		#downloadData(syn, sc3, testDataDir,replace)
		downloadData(syn, sc3, testDataDir)



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
