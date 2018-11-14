#! /usr/bin/env python
import dendropy
import sys
import os
from dendropy.calculate import treecompare 

def setuploop(treefilelist,outfile):

   largeDenseList = dendropy.TreeList();
   largeDenseTree = dendropy.Tree.get(path="largelengthDense.tt", schema="newick")
   largeDenseTreeString = largeDenseTree.as_string(schema="newick")
   largeDenseList.read(data=largeDenseTreeString, schema="newick")
   
   largeSparseList = dendropy.TreeList();
   largeSparseTree = dendropy.Tree.get(path="largelengthSparce.tt", schema="newick")
   largeSparseTreeString = largeSparseTree.as_string(schema="newick")
   largeSparseList.read(data=largeSparseTreeString, schema="newick")

   moderateDenseList = dendropy.TreeList();
   moderateDenseTree = dendropy.Tree.get(path="moderatelengthDense.tt", schema="newick")
   moderateDenseTreeString = moderateDenseTree.as_string(schema="newick")
   moderateDenseList.read(data=moderateDenseTreeString, schema="newick")

   moderateSparseList = dendropy.TreeList();
   moderateSparseTree = dendropy.Tree.get(path="moderatelengthSparce.tt", schema="newick")
   moderateSparseTreeString = moderateSparseTree.as_string(schema="newick")
   moderateDenseList.read(data=moderateSparseTreeString, schema="newick")

   smallDenseList = dendropy.TreeList();
   smallDenseTree = dendropy.Tree.get(path="smalllengthDense.tt", schema="newick")
   smallDenseTreeString = smallDenseTree.as_string(schema="newick')
   smallDenseList.read(data=smallDenseTreeString, schema="newick")

   smallSparseList = dendropy.TreeList();
   smallSparseTree = dendropy.Tree.get(path="smalllengthSparce.tt", schema="newick")
   smallSparseTreeString = smallSparseTree.as_string(schema="newick")
   smallSparseList.read(data=smallSparseTreeString, schema="newick")

   modelTrees = [largeDenseTree, largeSparseTree, moderateDenseTree, moderateSparseTree, smallDenseTree, smallSparseTree]
   NJTrees = []

   for file in enumerate(treefilelist):
      if "tre" not in file:
         continue
      NJTrees.append(file)

   for index, file in NJTrees

      if index < 20:
         treeList = dendropy.TreeList(taxon_namespace=largeDenseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(largeDenseTree, treeString)
         symdiff = symdiff/len(largeDenseTree.nodes())
         print(symdiff)

       if index < 40 and index > 19:
         treeList = dendropy.TreeList(taxon_namespace=largeSparseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(largeSparseTree, treeString)
         symdiff = symdiff/len(largeSparseTree.nodes())
         print(symdiff)

       if index < 60 and index > 39:
         treeList = dendropy.TreeList(taxon_namespace=moderateDenseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(moderateDenseTree, treeString)
         symdiff = symdiff/len(moderateDenseTree.nodes())
         print(symdiff)

       if index < 80 and index > 59:
         treeList = dendropy.TreeList(taxon_namespace=moderateSparseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(moderateSparseTree, treeString)
         symdiff = symdiff/len(moderateSparseTree.nodes())
         print(symdiff)

       if index < 100 and index > 79:
         treeList = dendropy.TreeList(taxon_namespace=smallDenseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(smallDenseTree, treeString)
         symdiff = symdiff/len(smallDenseTree.nodes())
         print(symdiff)

       if index < 120 and index > 99:
         treeList = dendropy.TreeList(taxon_namespace=smallSparseList.taxon_namespace);
         tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
         treeString = tree.as_string(schema="newick")
         treeList.read(data=treeString, schema="newick")
         symdiff = treecompare.weighted_robinson_foulds_distance(smallSparseTree, treeString)
         symdiff = symdiff/len(smallSparseTree.nodes())
         print(symdiff)

def main():
   directory = os.listdir('./NJtrees')
   setuploop(directory)
 
if __name__ == "__main__":
    main()
