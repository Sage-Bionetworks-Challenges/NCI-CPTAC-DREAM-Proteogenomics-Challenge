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
	if args.express:
		sc1 = "syn10902163"
		sc2 = "syn10902231"
		sc3 = "syn10902232"
	else:
		if args.round == '1':
			sc1 = "syn10807805"
			sc2 = "syn10617839"
			sc3 = 'syn10617842'
			replace = "_round_1"
		elif args.round == '2':
			sc1 = "syn10807805"
			sc2 = 'syn10617840'
			sc3 = 'syn10617843'
			replace = "_round_2"
		elif args.round == 'test':
			sc1 = "syn10807805"
			sc2 = 'syn10139559'
			sc3 = 'syn10139567'
		else:
			sc1 = "syn10807814"
			sc2 = 'syn10617841'
			sc3 = 'syn10617844'
			replace = "_final_round"

	if args.sc == 'sc1':
		downloadData(syn, sc1,testDataDir)
		removeFile = os.path.join(testDataDir, "data_true.txt")
		if args.express:
			for i in range(1,101):
				os.rename(os.path.join(testDataDir, "data_obs_%s.txt" % i), os.path.join(testDataDir, "data_test_obs_%s.txt" % i))
	elif args.sc == 'sc2':
		downloadData(syn, sc2, testDataDir)
		removeFile = os.path.join(testDataDir, "prospective_ova_pro_gold_express.txt")
	else:
		downloadData(syn, sc3, testDataDir)
		removeFile = os.path.join(testDataDir, "prospective_ova_phospho_gold_express.txt")
	
	if args.express:
		os.remove(removeFile)


if __name__ == '__main__':
	main()
