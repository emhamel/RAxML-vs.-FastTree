import os
import dendropy

f = open("/projects/tallis/hamel/RAxML-vs.-FastTree/InputData/results.txt", "w+")

taxa = dendropy.TaxonNamespace()
originalTrees = []
njTrees = []

conditions = ["smalllengthDense", "smalllengthSparse", "moderatelengthDense", "moderatelengthSparse", "largelengthDense", "largelengthSparse"]

originalTreesPath = "/projects/tallis/hamel/RAxML-vs.-FastTree/ModelTrees/"
njTreesPath = "/projects/tallis/hamel/RAxML-vs.-FastTree/ModelTrees/NJtrees/"

matching_files = [f for f in os.listdir(originalTreesPath) if ".tt" in f]
for file in matching_files:
	originalTrees.append(file)

for files in os.listdir(njTreesPath):
	njTrees.append(files)

scores = []
for i in range(0, len(conditions)):
	matching_files = [f for f in njTrees if conditions[i] in f]
	tree1 = dendropy.Tree.get(file=open(originalTreesPath + conditions[i] + ".tt"),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
	tree1.encode_bipartitions()
	f.write(conditions[i] + "/n")
	for tree in matching_files:
		tree2 = dendropy.Tree.get(file=open(njTreesPath + tree),
                	schema="newick", rooting='force-unrooted', taxon_namespace=taxa)
		tree2.encode_bipartitions()
		results = compare_trees.compare_trees(tree1, tree2)
		scores.add(results[5])
		f.write(str(results[5]) + "/n")
	
f.write("avr: " + str(sum(scores)/(len(scores))) + "/n")
