import os
import dendropy

f = open("/Users/emmahamel/Research/RAxMLvsFastTree/RF_results/results.txt", "w+")

taxa = dendropy.TaxonNamespace()
conditions = ["smallLengthDense", "smallLengthSparce", "moderateLengthDense", "moderateLengthSparce", "largeLengthDense", "largeLengthSparce"]

originalTreesPath = "/projects/tallis/hamel/RAxML-vs.-FastTree/ModelTrees/"
njTreesPath = "/projects/tallis/hamel/RAxML-vs.-FastTree/ModelTrees/NJtrees"

matching_files = [f for f in os.listdir(originalTreesPath) if ".tt" in f]
for file in matching_files:
	originalTrees.append(file)

for files in os.listdir(njTreesPath):
	njTrees.append(files)

for i in range(0, len(conditions)):
	matching_files = [f for f in njTrees if conditions[i] in f]
	tree1 = dendropy.Tree.get(file=open(originalTreesPath + "/" + conditions[i] + ".tt"),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
	tree1.encode_bipartitions()

	for tree in matching_files:
		tree2 = dendropy.Tree.get(file=open(njTreesPath + "/" + tree),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
		tree2.encode_bipartitions()
		results = compare_trees.compare_trees(tree1, tree2)
		print(results[5])
		f.write(str(results[5]) + "/n")
