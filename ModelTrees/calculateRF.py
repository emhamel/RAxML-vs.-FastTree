#! /usr/bin/env python
import dendropy
import sys
import os
from dendropy.calculate import treecompare 

def setuploop(treefilelist,outfile):

   largeDenseTree = dendropy.Tree.get(path="largelengthDense.tt", schema="newick")
   largeSparseTree = dendropy.Tree.get(path="largelengthSparce.tt", schema="newick")
   moderateDenseTree = dendropy.Tree.get(path="moderatelengthDense.tt", schema="newick")
   moderateSparseTree = dendropy.Tree.get(path="moderatelengthSparce.tt", schema="newick")  
   smallDenseTree = dendropy.Tree.get(path="smalllengthDense.tt", schema="newick")
   smallSparseTree = dendropy.Tree.get(path="smalllengthSparce.tt", schema="newick")

   modelTrees = [largeDenseTree, largeSparseTree, moderateDenseTree, moderateSparseTree, smallDenseTree, smallSparseTree]
   NJTrees = []

   for file in treefilelist:
      if "tre" not in file:
         continue
      tree = dendropy.Tree.get(path="./NJtrees/"+file, schema="newick")
      NJTrees.append(tree)

   for tree in NJTrees[0:19]:
      #print(largeDenseTree.as_string(schema="newick"),)
      #print(tree.as_string(schema="newick"),)
      symdiff = treecompare.symmetric_difference(largeDenseTree, tree)
      #symdiff = treecompare.weighted_robinson_foulds_distance(largeDenseTree, tree)
      #symdiff = largeDenseTree.robinson_foulds_distance(tree) / len(largeDenseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

   for tree in NJTrees[20:39]:
      symdiff = treecompare.weighted_robinson_foulds_distance(largeSparseTree, tree)
      #symdiff = largeSparseTree.robinson_foulds_distance(tree) / len(largeSparseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

   for tree in NJTrees[40:59]:
      symdiff = treecompare.weighted_robinson_foulds_distance(moderateDenseTree, tree)
      #symdiff = moderateDenseTree.robinson_foulds_distance(tree) / len(moderateDenseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

   for tree in NJTrees[60:79]:
      symdiff = treecompare.weighted_robinson_foulds_distance(moderateSparseTree, tree)
      #symdiff = moderateSparseTree.robinson_foulds_distance(tree) / len(moderateSparseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

   for tree in NJTrees[80:99]:
      symdiff = treecompare.weighted_robinson_foulds_distance(smallDenseTree, tree)
      #symdiff = smallDenseTree.robinson_foulds_distance(tree) / len(smallDenseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

   for tree in NJTrees[100:119]:
      symdiff = treecompare.weighted_robinson_foulds_distance(smallSparseTree, tree)
      #symdiff = smallSparseTree.robinson_foulds_distance(tree) / len(smallSparseTree.nodes())
      towrite2 = str(symdiff) + '\n'
      with open(outfile, 'a') as out2:
         out2.write(towrite2)

def main():
   directory = os.listdir('./NJtrees')
   outputFile = open('./RFdistances.txt', 'w')
   setuploop(directory,outputFile)
   file.close()         
 
if __name__ == "__main__":
    main()
