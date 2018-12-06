import dendropy
import compare_trees
import os

originalTrees = []
njTrees = []
taxa = dendropy.TaxonNamespace()

f = open("/Users/emmahamel/Research/RAxMLvsFastTree/InputData/results_round3.txt", "w+")

conditions = ["small_dense", "small_sparse", "moderate_dense", "moderate_sparse", "large_dense", "large_sparse"]

originalTreesPath = "/Users/emmahamel/Research/RAxMLvsFastTree/ModelTrees"
njTreesPath = "/Users/emmahamel/Research/RAxMLvsFastTree/ModelTrees/NJtrees"

matching_files = [f for f in os.listdir(originalTreesPath) if ".tt" in f]
for file in matching_files:
	originalTrees.append(file)


for files in os.listdir(njTreesPath):
	njTrees.append(files)

scores = 0
for i in range(0, len(conditions)):
	matching_files = [f for f in njTrees if conditions[i] in f]
	tree1_path = originalTreesPath + "/" + conditions[i] + ".tt"
	tree1 = dendropy.Tree.get(file=open(tree1_path),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
	tree1.encode_bipartitions()

	for tree in matching_files:
		tree2_path = njTreesPath + "/" + tree
		tree2 = dendropy.Tree.get(file=open(tree2_path),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
		tree2.encode_bipartitions()
		results = compare_trees.compare_trees(tree1, tree2)
		f.write(tree1_path + "," + tree2_path + "," + conditions[i] + "," + str(results[5]) + "\n")
		scores += results[5]

	print(scores)
	f.write("avr: " + str(scores/20) + "\n")
	scores = 0
