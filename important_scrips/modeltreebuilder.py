import dendropy
import random

treeDictionaryDenseTaxa = {}
treeDictionarySparceTaxa = {}

work='/projects/tallis/hamel/RAxML-vs.-FastTree/ModelTrees/'

for i in range (0,3):
	treeDictionaryDenseTaxa[i] = dendropy.Tree.get(file=open(work + 'Originaltree.tt'),
                	schema="newick")
	treeDictionarySparceTaxa[i] = dendropy.Tree.get(file=open(work + 'Originaltree.tt'),
                	schema="newick")

setConstant = [10, 50, 100]

for k in range(0, 3):
	treeDictionaryDenseTaxa[k].scale_edges(setConstant[k])
	treeDictionarySparceTaxa[k].scale_edges(setConstant[k])

randomNodes = []
nodesArray = treeDictionarySparceTaxa[0].leaf_nodes()

for p in range (0, 50):
	choice = random.choice(nodesArray)
	randomNodes.append(choice.taxon.label)
	nodesArray.remove(choice)

for l in range(0, len(treeDictionarySparceTaxa)):
	treeDictionarySparceTaxa[l] = treeDictionarySparceTaxa[l].extract_tree_with_taxa_labels(randomNodes)

setName = ["small", "moderate", "large"]

for i in range (0,3):
	treeDictionaryDenseTaxa[i].write(path=work + setName[i] + '_dense.tt', 
		schema="newick")
	treeDictionarySparceTaxa[i].write(path=work + setName[i] + '_sparse.tt', 
		schema="newick")
