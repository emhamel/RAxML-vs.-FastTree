#! /usr/bin/env python
import dendropy
import sys
import os
from dendropy.calculate import treecompare 

def setuploop(treefilelist):

   tns = dendropy.TaxonNamespace()

   largeDenseTree = dendropy.Tree.get_from_path("largelengthDense.tt", "newick", taxon_namespace=tns)
   largeSparseTree = dendropy.Tree.get_from_path("largelengthSparce.tt", "newick", taxon_namespace=tns)
   moderateDenseTree = dendropy.Tree.get_from_path("moderatelengthDense.tt", "newick", taxon_namespace=tns)
   moderateSparseTree = dendropy.Tree.get_from_path("moderatelengthSparce.tt", "newick", taxon_namespace=tns)
   smallDenseTree = dendropy.Tree.get_from_path("smalllengthDense.tt", "newick", taxon_namespace=tns)
   smallSparseTree = dendropy.Tree.get_from_path("smalllengthSparce.tt", "newick", taxon_namespace=tns)

   largeDenseTree.encode_bipartitions()
   largeSparseTree.encode_bipartitions()
   moderateDenseTree.encode_bipartitions()
   moderateSparseTree.encode_bipartitions()
   smallDenseTree.encode_bipartitions()
   smallSparseTree.encode_bipartitions()

   NJTrees = [x for x in treefilelist if "tre" in x]
  
   for file in NJTrees:
      if NJTrees.index(file) < 20:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(largeDenseTree, tree2)
         rf_error = rf/2*len(largeDenseTree.internal_edges())
         print(rf_error)

      if NJTrees.index(file) < 40 and NJTrees.index(file) > 19:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(largeSparseTree, tree2)
         rf_error = rf/2*len(largeSparseTree.internal_edges())
         print(rf_error)

      if NJTrees.index(file) < 60 and NJTrees.index(file) > 39:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(moderateDenseTree, tree2)
         rf_error = rf/2*len(moderateDenseTree.internal_edges())
         print(rf_error)

      if NJTrees.index(file) < 80 and NJTrees.index(file) > 59:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(moderateSparseTree, tree2)
         rf_error = rf/2*len(moderateSparceTree.internal_edges())
         print(rf_error)

      if NJTrees.index(file) < 100 and NJTrees.index(file) > 79:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(smallDenseTree, tree2)
         rf_error = rf/2*len(smallDenseTree.internal_edges())
         print(rf_error)

      if NJTrees.index(file) < 120 and NJTrees.index(file) > 99:
         tree2 = dendropy.Tree.get_from_path("./NJTrees/"+file, "newick", taxon_namespace=tns)
         tree2.encode_bipartitions()
         rf = treecompare.weighted_robinson_foulds_distance(smallSparseTree, tree2)
         rf_error = rf/2*len(smallSparceTree.internal_edges())
         print(rf_error)

def main():
   directory = os.listdir('./NJtrees')
   setuploop(directory)
 
if __name__ == "__main__":
    main()
