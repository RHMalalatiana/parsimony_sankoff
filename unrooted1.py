#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 13:33:15 2022

@author: hanitrinialamr
"""

def generate_unrooted_trees(num_of_species):
  n=num_of_species
  remaining_species=list(range(2,n+1))
  lst_out=[]
  for subtree in generate_subtrees(remaining_species):
    new_tree=[1,subtree]
    lst_out.append(new_tree) 
  return(lst_out)

def generate_subtrees(remaining_species):
  if len(remaining_species) == 1:
    return([remaining_species[0]])
  else : 
    first=min(remaining_species) 
    remaining_species.remove(first)  # check this
    out_list=[] 
    for left,right in partitions(remaining_species):
      if len(right) != 0:
        left=[first]+left
        for left_subtree in generate_subtrees(left): 
          for right_subtree in generate_subtrees(right):
             out_list.append([left_subtree,right_subtree])
    return(out_list) 

def partitions(lst):
  length=len(lst)
  partitions_out=[]
  num_partitions=2**length
  for part_int in range(num_partitions):
    left =[]
    right=[]
    for k in range(length):
      modulus=2**(k+1)
      if part_int % modulus == 0 :
        left.append( lst[k])
      else : 
        right.append( lst[k] )
        part_int=part_int//2
    new_partition=(left,right)
    partitions_out.append(new_partition)
  return(partitions_out)

#r=generate_subtrees([1,2,3])  
#print(r) 
#
#r=generate_unrooted_trees(5)  
#print(r) 

def replace_tree_labels(tree,label_list):
  left,right=tree
  if type(left) == list :
    left=replace_tree_labels(left,label_list)
  else :
    left=label_list[left-1]
  if type(right) == list :
    right=replace_tree_labels(right,label_list)
  else :
    right=label_list[right-1]
  return([left,right])

def replace_tree_labels_in_tree_list(tree_list,label_list):
  output=[]
  for tree in tree_list:
    output.append(
       replace_tree_labels(tree,label_list)
    )
  return(output)  

def unroot_generate_driver(label_list):
  tree_list= generate_unrooted_trees(len(label_list))  
  return  (replace_tree_labels_in_tree_list(tree_list,label_list))

#label_list=[ 'a', 'b', 'c', 'd', 'e']
#rr=replace_tree_labels_in_tree_list(r,label_list) 
#print(rr) 
#
#label_list=[ 'a', 'b', 'c', 'd', 'e']
#rrr=unroot_generate_driver(label_list)
#print(rrr)


