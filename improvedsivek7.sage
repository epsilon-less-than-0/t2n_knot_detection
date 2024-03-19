### Code for  "Khovanov homology and the cinquefoil", arXiv:2105.12102
###
### Enumerate Stallings 7-braids for which the braid axis might
### lift in the branched double cover to a knot with the same
### Khovanov homology as T(2,7).
###
### The function check_stallings_braids() does the following:
### - Find all Stallings 7-braids b with the following property:
###     b is quasipositive, and in the branched double cover of
###     the closure of b, the braid axis lifts to a knot A with
###     the same Alexander polynomial as T(2,7).
### - Verify that for all such b, either the lift A is T(2,7),
###     or b is not exchangeable because the closure of b^4
###     is not fibered (its Alexander polynomial is not monic).
import time
start_time = time.time()
def elembraids(n=7,positive=True):
  ## Return a tuple of the elementary braids on n strands.
  BG = BraidGroup(n)
  eblist = []
  for i in range(1,n):
    for j in range(i+1,n+1):
      conjugator = BG( list(range(j-1,i,-1)) )
      s = conjugator * BG([i]) * conjugator**(-1) ## exchange strands i,j
      if positive:
        eblist.append(s) ## positive generators only
      else:
        eblist += [s, s**(-1)]
  return tuple(eblist)

def stallings(n=7,positive=True,verbose=False):
  eblist = elembraids(n,positive)
  answers = []
  itercount = 0
  for braidtuple in cartesian_product_iterator( [eblist]*(n-1) ):
    itercount += 1
    b = prod(braidtuple)
    if b.components_in_closure() == 1: ## it's a knot, so it's an unknot
      A = b.burau_matrix(reduced=True)
      ## chech Alexander poly same as T(2,7) and also braid is pseudoanosov
      if bool(A.subs(t=-1).characteristic_polynomial() == x^6 - x^5 + x^4 - x^3 + x^2 - x + 1) == True and b.is_pseudoanosov() == True: 
        answers.append(b)
    if verbose and itercount%1000 == 0:
      print("%d braids checked out of %d, %d successful"%(itercount,len(eblist)**(n-1),len(answers)))
      print("--- %s seconds ---" % (time.time() - start_time))

  return answers

def check_stallings_braids(verbose=True):
  ### Prove that if b is a pseudo-Anosov, exchangeable
  ### 7-braid, then in the branched double cover of the
  ### closure of b, the braid axis lifts to the knot
  ### fakeT27.  Then compute its knot Floer homology
#   def check_isometric(M1,M2):
    ### SnapPy sometimes randomly fails to find isometries,
    ### so we repeat until it (hopefully) does
    # for _ in range(10):
    #   try:
    #     if M1.is_isometric_to(M2):
    #       return True
    #   except RuntimeError:
    #     M1.randomize()
    #     M2.randomize()
    # return False

  print("Generating the Stallings 7-braids with correct Alexander polynomial and are pseudoanosov...")
  sb_list_all = stallings(7,True,verbose)
  print("%d candidate braids found"%(len(sb_list_all)))
  print("--- %s seconds ---" % (time.time() - start_time))


  # Take one braid representative of each conjugacy class
  sb_list = [BraidGroup(7)([6,5,4,3,2,1])] # start with a braid giving T(2,7)
  itercount = 1
  for b in sb_list_all:
    itercount += 1
    if not True in [b.is_conjugated(beta) for beta in sb_list]:
      sb_list.append(b)
    print("%d candidate braids checked out of %d, %d added to distinct conjugacy classes"%(itercount,len(sb_list_all),len(sb_list)))
  print("%d distinct conjugacy classes found, including %s"%(len(sb_list),str(sb_list[0])))
  print("--- %s seconds ---" % (time.time() - start_time))

  print("Checking that 4th powers of braids are not fibered...")
  still_successful = True
  for sb in sb_list[1:]: # check all the classes that don't lift to T(2,7) because the 0th index is T(2,7)
    ap = (sb**4).alexander_polynomial()
    if abs(ap.coefficients()[0]) == 1:
      print("*** braid might be exchangeable:")
      print("***    b=%s"%str(sb.Tietze()))
      still_successful = False

  if still_successful == False:
    print("--- %s seconds ---" % (time.time() - start_time))
    return False
  else: 
    print("All exchangeable braid axes lifted to T(2,7).")
    print("--- %s seconds ---" % (time.time() - start_time))
    return True    
  
print(check_stallings_braids())