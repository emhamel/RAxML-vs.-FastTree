import dendropy
import random

treeDictionaryDenseTaxa = {}
treeDictionarySparseTaxa = {}

for i in range (0,3):
	treeDictionaryDenseTaxa[i] = dendropy.Tree.get(file=open('/users/emmahamel/Research/ModelTrees/Originaltree.tt'),
                	schema="newick")
	treeDictionarySparseTaxa[i] = dendropy.Tree.get(file=open('/users/emmahamel/Research/ModelTrees/Originaltree.tt'),
                	schema="newick")

setConstant = [.2, 1, 5]

for k in range(0, 3):
	treeDictionaryDenseTaxa[k].scale_edges(setConstant[k])
	treeDictionarySparseTaxa[k].scale_edges(setConstant[k])

randomNodes = []
nodesArray = treeDictionarySparseTaxa[0].leaf_nodes()

for p in range (0, 50):
	choice = random.choice(nodesArray)
	randomNodes.append(choice.taxon.label)
	nodesArray.remove(choice)

for l in range(0, len(treeDictionarySparseTaxa)):
	treeDictionarySparceTaxa[l] = treeDictionarySparceTaxa[l].extract_tree_with_taxa_labels(randomNodes)

setName = ["small", "moderate", "large"]

for i in range (0,3):
	treeDictionaryDenseTaxa[i].write(path='/users/emmahamel/Research/ModelTrees/' + setName[i] + 'lengthDense.tt', 
		schema="newick")
	treeDictionarySparceTaxa[i].write(path='/users/emmahamel/Research/ModelTrees/' + setName[i] + 'lengthSparce.tt', 
		schema="newick")
