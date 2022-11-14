
import copy 

def get_leaf_list(tree):
  if is_leaf(tree):
    return([tree])
  left,right=tree
  return( get_leaf_list(left)
         +get_leaf_list(right)
  )

def val(tree):
  leaf_list=get_leaf_list(tree)
  #print(leaf_list)
  return( min(leaf_list))

def is_leaf(tree):
  return( type(tree) != list )

def reorder(tree):
   left,right=tree
   if val(left) > val(right):
     left,right=right,left
   if not is_leaf(left):
     left=reorder(left)
   if not is_leaf(right):
     right=reorder(right)
   return([left,right])

def reroot(tree):
   left,right=tree
   if val(left) > val(right):
     left,right=right,left
   if is_leaf(left):
     return([left,right])
   else :
     ll,lr=left
     if val(ll) > val(lr):
       ll,lr=lr,ll
     tree=[ll,[lr,right]]
     return(reroot(tree)) 

def canonical(tree):
   tree=reroot(tree) 
   tree=reorder(tree)
   return(tree) 

def canonical_rooted_list(tree_list):
  ctree_list=[ canonical(tree) for tree in tree_list]
  result=remove_duplicates(ctree_list)
  return(result)

def remove_duplicates(lst_in):
  lst=copy.deepcopy(lst_in)
  lst.reverse()
  lst_out=[lst.pop()]
  try:
    while(True):
      ele=lst.pop()
      is_duplicate=False
      for ele1 in lst_out:
        if ele == ele1 :
          is_duplicate=True
      if is_duplicate==False :
        lst_out.append(ele) 
  except IndexError:
    pass
  return(lst_out) 

        
  
