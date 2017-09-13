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
	parser.add_argument("round", help="Round", choices=['1','2',"final",'test'])
	parser.add_argument("--express", help="Express lane", action="store_true", default=False)
	args = parser.parse_args()
	syn = synapseclient.login()
	downloadDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"goldstandard")
	if not os.path.exists(downloadDir):
		os.mkdir(downloadDir)
	
	# downloadData(syn, "syn10165897",downloadDir)
	# downloadData(syn, "syn10165898",downloadDir)
	# downloadData(syn, "syn10165900",downloadDir)
	downloadData(syn, "syn10164401",downloadDir)
	sc1 = syn.get("syn10164502")
	if args.round == '1':
		sc2 = syn.get("syn10763208")
		sc3 = syn.get("syn10763237")
	elif args.round == '2':
		sc2 = syn.get("syn10763217")
		sc3 = syn.get("syn10763243")
	elif args.round == 'test':
		sc2 = syn.get("syn10514976")
		sc3 = syn.get("syn10666694")
	else:
		sc2 = syn.get("syn10763225")
		sc3 = syn.get("syn10763252")
	shutil.copy(sc1.path, downloadDir)
	shutil.copy(sc2.path, os.path.join(downloadDir, "prospective_ova_pro_gold.txt"))
	shutil.copy(sc3.path, os.path.join(downloadDir, "prospective_ova_phospho_gold.txt"))
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
