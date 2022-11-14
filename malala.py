import copy
import random
import math as M
import numpy as np 
from canonical import *

#  TREE GENERATOR FUNCTIONS

""""
Note that both the functions below are recursive.
This is a useful programming trick when dealing with trees.
"""

def trees_with_extra_species(input_tree,new_species):
    '''
    Inputs: 
      input_tree - in the form of a list of lists (see above)
      new_species - N.B. can be any type other than a list

    Output: list of all possible trees resulting from adding the new_species 
            to each branch of input_tree
    '''
    tree_list_out=[]
    if type(input_tree[0])!=list and type(input_tree[1])==list:
        tree_list_out.append([[input_tree[0],new_species],input_tree[1]])
        tree_list_out.append([[input_tree[0],input_tree[1]],new_species])
        for i in trees_with_extra_species(input_tree[1],new_species):
            tree_list_out.append([input_tree[0],i])
    elif type(input_tree[1])!=list and type(input_tree[0])==list:
        tree_list_out.append([input_tree[0],[input_tree[1],new_species]])
        tree_list_out.append([[input_tree[0],input_tree[1]],new_species])
        for i in trees_with_extra_species(input_tree[0],new_species):
            tree_list_out.append([i,input_tree[1]])
    elif type(input_tree[0])!=list and type(input_tree[1])!=list:
        tree_list_out.append([[input_tree[0],new_species],input_tree[1]])
        tree_list_out.append([input_tree[0],[input_tree[1],new_species]])
        tree_list_out.append([[input_tree[0],input_tree[1]],new_species])
    else:
        tree_list_out.append([input_tree,new_species])
        for new_subtree in trees_with_extra_species(input_tree[0],new_species):
            tree_list_out.append([new_subtree,input_tree[1]])
        for new_subtree in trees_with_extra_species(input_tree[1],new_species):
            tree_list_out.append([input_tree[0],new_subtree])     
    return(tree_list_out)  

def enumerate_labelled_trees(leaf_list):
    '''
    Input: a list of leaf labels (typically enumerated using integers)
    Output: a list of all labelled bifurcating trees using all labels in input list
    '''
    all_trees=[[leaf_list[0],leaf_list[1]]]
    for new_leaf in leaf_list[2:]:
      new_all_trees=[]
      for tree in all_trees:
        new_all_trees+=trees_with_extra_species(tree,new_leaf)
      all_trees=new_all_trees
    return(all_trees)



# COUNTING CHANGE (SANKOFF)

def evaluate_Sankoff_tableau(sankoff_tableau_left,sankoff_tableau_right, cost_matrix):
    '''
    Input: Sankoff tableaux for the left and right children of node
    Output: Sankoff tableaux of node 
    '''
    dm=len(sankoff_tableau_left)
    sTableauOut=np.zeros(dm)
    l_vec  =np.zeros(dm)
    r_vec  =np.zeros(dm)
    for i in range(dm):
        l_vec=sankoff_tableau_left +cost_matrix[i,:]
        r_vec=sankoff_tableau_right+cost_matrix[i,:]
        sTableauOut[i]=min(l_vec)+min(r_vec)
    return(sTableauOut)

def Sankoff(tree, alphabet, speciesCharacters, cost_matrix):
    '''
    This function computes the number of changes for one character per species
    Input:
        - tree
        - alphabet
        - speciesCharacters
        - cost_matrix
    Output:
        - score, vector cost of the root
    '''
    def leaf_tableau(char_in):
      index=alphabet.index(char_in)
      tableau = [float('inf')] * len(alphabet)
      tableau[index]=0. 
      return(tableau)
    score=0
    if type(tree[0])!=list and type(tree[1])==list:
        score,root_tableau=Sankoff(tree[1],alphabet,speciesCharacters,cost_matrix)
        lch=speciesCharacters[tree[0]-1]
        l_tableau=leaf_tableau(lch)
        root_tableau= evaluate_Sankoff_tableau(l_tableau,root_tableau,cost_matrix)
    elif type(tree[1])!=list and type(tree[0])==list:
        score,root_tableau=Sankoff(tree[0],alphabet,speciesCharacters,cost_matrix)
        rch=speciesCharacters[tree[1]-1]
        r_tableau=leaf_tableau(rch)
        root_tableau= evaluate_Sankoff_tableau(root_tableau,r_tableau,cost_matrix)
    elif type(tree[0])==list and type(tree[1])==list:
        score,root_tableau_l=Sankoff(tree[0],alphabet,speciesCharacters,cost_matrix)
        score,root_tableau_r=Sankoff(tree[1],alphabet,speciesCharacters,cost_matrix)
        root_tableau= evaluate_Sankoff_tableau(root_tableau_l,root_tableau_r,cost_matrix)
    else:    #    both children are leaves 
        lch=speciesCharacters[tree[0]-1]
        l_tableau=leaf_tableau(lch)
        rch=speciesCharacters[tree[1]-1]
        r_tableau=leaf_tableau(rch)
        root_tableau= evaluate_Sankoff_tableau(l_tableau,r_tableau,cost_matrix)
    score=min(root_tableau)
    return((score, root_tableau))   

def Sankoff_tree(tree, data, alphabet, cost_matrix):
    '''
    This function counts the number of changes for all characters of species
    Input:
            - tree
            - data : matrix which represents all the characters of each species, species in the row and character in the column    
                 (actually a list of species genomes, where the species genome is a list of characters)
            - alphabet 
            - cost_matrix
    Output:
            - count: parsimony score of input tree
    '''
    count=0
    nbases=len(data[0])
    for ibase in range(nbases):
        baseList=[species[ibase] for species in data]
        score,root_tableau=Sankoff(tree,alphabet,baseList,cost_matrix)
        count+=score
    return(count)


# MOST PARSIMONIOUS TREE (SANKOFF)

def parsimonious_Sank(tree_list,data,alphabet,cost_matrix):
    '''
    Input: 
            - tree_list = all possible trees for a given set of species (or a subset thereof)
            - data: matrix which represent all the Character of each species, species in the row and character in the column    
            - alphabet : alphabet  of allowed character states
            - cost: cost matrix
    Output: most parsimonious unrooted tree
    '''
    minim=Sankoff_tree(tree_list[0],data,alphabet,cost_matrix) 
    most_parsimonious=[tree_list[0]]
    for tree in tree_list[1:]:
        score=Sankoff_tree(tree,data,alphabet,cost_matrix)
        if score<minim:
            most_parsimonious=[tree]
            minim=score
        elif score==minim:
            most_parsimonious.append(tree)
    return((canonical_rooted_list(most_parsimonious[0]), minim))


# SIMULATED DATA GENERATION

def jukesCantor(genome,gen_time):
  """
  Input: genome, a character string consisting of a, c, g, and t
         gen_time, genetic time
  Output: a randomly modified string according to the Juke's cantor model
    The following coupled ODEs are integrated from t=0 to t=gen_time

    dot a = -a + (c + g + t) /3
    dot c = -a + (g + t + a) /3
    dot g = -g + (t + a + c) /3
    dot t = -t + (a + c + g) /3
  
    using the solution 
 
    a(t)=(3/4)*exp(-t)*(a(0)-(c(0)+g(0)+t(0))/3)
        +(1/4)*(a(0)+c(0)+g(0)+t(0))

    and similarly for c(t), g(t), and t(t).

    In other words, if 'a' is the initial base, its persistence
    probability after a time t is 

    3/4 exp(-t) + 1/4

    and likewise for the other three base inputs.

  """
  if type(gen_time) !=  float :
    raise Exception("Time must be a float")
  for b in genome:
    if not ( b=='a' or b=='c' or b=='g' or b=='t' ):
       raise Exception("Invalid input genome")
  persistence_prob=(3.*M.exp(-gen_time)+1.)/4.
  new_genome=[]
  for base in genome:
    p=random.random()
    if p < persistence_prob:
      new_base=base
    else: 
      if base=='a':
        seq=('c','g','t')
      if base=='c':
        seq=('a','g','t')
      if base=='g':
        seq=('a','c','t')
      if base=='t':
        seq=('a','c','g')
      new_base=random.choice(seq)
    new_genome.append(new_base)
  return(new_genome)

# The following routine, which applies mutations randomly, chooses the Jukes-Cantor algorithm,
# but some other model can be used by changing the function below.

def evolve(genome,time):
   return(jukesCantor(genome,time))

def generateHelper(initialGenome,listIn):
 time=listIn[0]
 newGenome=evolve(initialGenome,time)
 listIn[0]=newGenome
 if not ( len(listIn) == 1 or len(listIn) == 3):
    raise Exception("Lists must have one or three elements") 
 if len(listIn) == 1:
   return
 else:
   generateHelper(newGenome,listIn[1])
   generateHelper(newGenome,listIn[2])
   return 

def generateDriver(initialGenome,listIn):
  """
  This algorithm assumes a tree with input of the form:

  [ genTime, leftBranch, rightBranch]

  where the branches can either be tripleton lists of the format above
  or leaves represtented as a single nonnegative number enclosed in a list.
  genTime above is a non-negative real number. A deep copy of the input list is made,
  and the genetic times are replaced with the stochastically evolved genomes.

  An initial genone at the root of the binary tree is required as input, and this
  must be in a format of a list of base charcters, i.e., 'a', 'c', 'g', or 't'.

  """
  listInBis=copy.deepcopy(listIn)
  generateHelper(initialGenome,listInBis)
  return(listInBis)
   
