import sys, os, csv
from ete3 import Tree

def main(args):
	labels = {}
	# load tab
	with open(args.labels_tab, "r") as f2:
		reader=csv.reader(f2,delimiter='\t')
		for row in reader:
			if row[0] != "id":
				if args.use_labels == "label2":
					labels[row[0]] = row[4]
				elif args.use_labels == "label3":
					labels[row[0]] = row[5]
				elif args.use_labels == "label1":
					labels[row[0]] = row[3]
				elif args.use_labels == "organism":
					labels[row[0]] = row[2]
				else: # ncbi_taxid
					taxid = row[1].split(" ")[0]
					labels[row[0]] = taxid

	keys = labels.keys()

	#  load tree
	t = Tree(args.input_tree,format=1, quoted_node_names=False)

	for node in t.traverse():
		if node.is_leaf() is True and node.name in keys:
			node.name = labels[node.name]

	# write tree to a file
	t.write(format=1, outfile="{output_file}".format(output_file=args.output_file))

if __name__ == "__main__" and '__file__' in globals():
	import argparse
	cpu_count=os.cpu_count()
	parser = argparse.ArgumentParser(prog="label",description="label sequences")
	parser.add_argument('-i','--input_tree', dest="input_tree", required=True, help='input tree file')
	parser.add_argument('-l','--labels_tab', dest="labels_tab", required=True, help='tabular file with labels')
	parser.add_argument('--use_labels', dest="use_labels", required=False, type=str.lower,choices=["label1", "label2", "label3","organism","ncbi_taxid"], default="label3", help='column to use')
	parser.add_argument('-o','--output_file', dest="output_file", required=True, help="output folder and base filename")

	if len(sys.argv) == 1:
		sys.stderr.write("No arguments provided, printing help menu ...\n")
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	main(args)
