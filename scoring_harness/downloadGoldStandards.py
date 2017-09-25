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
	parser.add_argument("round", help="Round", choices=['1','2',"final"])
	args = parser.parse_args()
	syn = synapseclient.login()
	downloadDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"goldstandard")
	expressDir = os.path.join(downloadDir,"express")
	if not os.path.exists(downloadDir):
		os.mkdir(downloadDir)
	if not os.path.exists(expressDir):
		os.mkdir(expressDir)
	#sc1 test, round1,2
	downloadData(syn, "syn10807805",downloadDir)
	downloadData(syn, "syn10807065",downloadDir)
	sc2_test = syn.get("syn10514976")
	sc3_test = syn.get("syn10666694")
	if args.round == '1':
		sc2 = syn.get("syn10763208")
		sc3 = syn.get("syn10763237")
	elif args.round == '2':
		sc2 = syn.get("syn10763217")
		sc3 = syn.get("syn10763243")
	else:
		downloadData(syn,"syn10807068",downloadDir)
		downloadData(syn,"syn10807814",downloadDir)
		sc2 = syn.get("syn10763225")
		sc3 = syn.get("syn10763252")

	shutil.copy(sc2.path, os.path.join(downloadDir, "prospective_ova_pro_gold.txt"))
	shutil.copy(sc3.path, os.path.join(downloadDir, "prospective_ova_phospho_gold.txt"))
	shutil.copy(sc2_test.path, os.path.join(downloadDir, "prospective_ova_pro_gold_complete.txt"))
	shutil.copy(sc3_test.path, os.path.join(downloadDir, "prospective_ova_phospho_gold_complete.txt"))

	#Express lane data
	downloadData(syn, "syn10902163",expressDir)
	for i in range(1,101):
		shutil.copy(os.path.join(expressDir,"data_true.txt"),os.path.join(expressDir,"data_test_true_%s.txt" % i)) 
	sc2 = syn.get("syn10903693")
	shutil.copy(sc2.path, os.path.join(expressDir, "prospective_ova_pro_gold_express.txt"))
	sc3 = syn.get("syn10903614")
	shutil.copy(sc3.path, os.path.join(expressDir, "prospective_ova_phospho_gold_express.txt"))

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
