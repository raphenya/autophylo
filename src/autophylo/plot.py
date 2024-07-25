from ete3 import Tree, TreeStyle, TextFace, NodeStyle, RectFace
import os, sys, csv

def set_color(node,color):
	node.add_feature("color",color)
	node.add_face(TextFace(node.name + "  ", fsize=80, ftype='Arial', tight_text=True), column=0)
	# set background color of node
	style = NodeStyle()
	style["bgcolor"] = color
	style["size"] = 10
	node.set_style(style)

def get_phyla_taxid(node_phyla, needle):
	for k in node_phyla:
		if needle in k:
			return k

def main(args):
	t = Tree(args.input_tree,format=0, quoted_node_names=True)
	# read node.name = phyla
	node_phyla = {}
	with open(args.labels_tab, "r") as f2:
		reader=csv.reader(f2,delimiter='\t')
		for row in reader:
			if row[0] != "id":
				node_phyla[row[5]] = row[6]

	for node in t.traverse("postorder"):
		if node.name != "":
			i = get_phyla_taxid(node_phyla, node.name)
			try:
				if args.putative_label in node.name: # putative
					set_color(node,"#6a0dad")#FF0000 #EE7733
				elif node_phyla[i] in "1224": # Proteobacteria
					set_color(node,"#CC3311")
				elif node_phyla[i] in "1239": # Firmicutes
					set_color(node,"#009988")
				elif node_phyla[i] in "201174": # Actinobacteria
					set_color(node,"#0077BB")
				elif node_phyla[i] in "203490": # Fusobacteriia
					set_color(node,"#FFEE99")
				elif node_phyla[i] in "32066": # Fusobacteria
					set_color(node,"#FFEE99")
				elif node_phyla[i] in "508458": # Synergistetes
					set_color(node,"#EE3377")
				elif node_phyla[i] in "74201": # Verrucomicrobia
					set_color(node,"#BBBBBB")
				elif node_phyla[i] in "976": # Bacteroidetes
					set_color(node,"#33BBEE")
				else: # Grey
					set_color(node,"#FFFFFF")
			except Exception as e:
				print("Exception: ", e, node.name, " => ", i)
		else:
			node.add_face(TextFace(node.support, fsize=80, ftype='Arial', tight_text=False), column=0)

	ts = TreeStyle()

	# Draw a tree
	ts.mode = "c"

	ts.scale =  50 #20 # 120 pixels per branch length unit
	ts.branch_vertical_margin = 1000 # 10 pixels between adjacent branches

	# We will add node names
	ts.show_leaf_name = False

	# Show branch data
	ts.show_branch_length = False
	ts.show_branch_support = False # True

	ts.branch_vertical_margin =  2
	# ts.draw_aligned_faces_as_table = True
	# ts.root_opening_factor=1

	ts.optimal_scale_level = 'full'
	ts.allow_face_overlap = False
	ts.arc_start = -180
	ts.arc_span = 360

	ts.title.add_face(TextFace("AutoPhylo plot", fsize=20), column=0)

	# Bacteroidetes
	ts.legend.add_face(RectFace(200,200,"black","#33BBEE"), column=0)
	ts.legend.add_face(TextFace("Bacteroidetes" + "  ", fsize=80, ftype='Arial', tight_text=True), column=1)

	# Proteobacteria
	ts.legend.add_face(RectFace(200,200,"black","#CC3311"), column=2)
	ts.legend.add_face(TextFace("Proteobacteria" + "  ", fsize=80, ftype='Arial', tight_text=True), column=3)

	# Firmicutes
	ts.legend.add_face(RectFace(200,200,"black","#009988"), column=4)
	ts.legend.add_face(TextFace("Firmicutes" + "  ", fsize=80, ftype='Arial', tight_text=True), column=5)

	# Actinobacteria
	ts.legend.add_face(RectFace(200,200,"black","#0077BB"), column=6)
	ts.legend.add_face(TextFace("Actinobacteria" + "  ", fsize=80, ftype='Arial', tight_text=True), column=7)

	# Fusobacteria
	ts.legend.add_face(RectFace(200,200,"black","#FFEE99"), column=0)
	ts.legend.add_face(TextFace("Fusobacteria" + "  ", fsize=80, ftype='Arial', tight_text=True), column=1)

	# Synergistetes
	ts.legend.add_face(RectFace(200,200,"black","#EE3377"), column=2)
	ts.legend.add_face(TextFace("Synergistetes" + "  ", fsize=80, ftype='Arial', tight_text=True), column=3)

	# Verrucomicrobia
	ts.legend.add_face(RectFace(200,200,"black","#BBBBBB"), column=4)
	ts.legend.add_face(TextFace("Verrucomicrobia" + "  ", fsize=80, ftype='Arial', tight_text=True), column=5)

	# Putative
	ts.legend.add_face(RectFace(200,200,"black","#6a0dad"), column=6)
	ts.legend.add_face(TextFace("Putative" + "  ", fsize=80, ftype='Arial', tight_text=True), column=7)

	p = os.path.join("{}.svg".format(args.output_file))
	# write tree
	t.render(p, tree_style=ts, w=1000, h=1080, dpi=600)

if __name__ == "__main__" and '__file__' in globals():
	import argparse
	cpu_count=os.cpu_count()
	parser = argparse.ArgumentParser(prog="plot-tree",description="plot tree")
	parser.add_argument('-i','--input_tree', dest="input_tree", required=True, help='input tree file')
	parser.add_argument('--putative_label', dest="putative_label", required=False, help='specify putative label to highlight')
	parser.add_argument('-l','--labels_tab', dest="labels_tab", required=True, help='tabular file with labels')
	parser.add_argument('-o','--output_file', dest="output_file", required=True, help="output folder and base filename")

	if len(sys.argv) == 1:
		sys.stderr.write("No arguments provided, printing help menu ...\n")
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()

	main(args)
