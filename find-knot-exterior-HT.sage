#Run through the Snappy CensusKnots database of knot exteriors
from snappy import *
b=BraidGroup(7)([1, 2, 1, -2, 3, 2, 1, -2, -3, 4, 3, 2, 1, -2, -3, -4, 5, 4, 3, 2, -3, -4, -5, 6, 5, 4, 3, -4, -5, -6])
M = Manifold('Braid%s(2,0)(0,0)'%(list(b.Tietze())))
homologyZZ = lambda Y: Y.homology().elementary_divisors() == [0]
coverlist = list(filter(homologyZZ, M.covers(2)))
M = coverlist[0]

def check_knot_exterior():
    manifoldfound = False
    exterior_list = []
    for N in HTLinkExteriors:
        if M.is_isometric_to(N) == True:
            exterior_list.append(N)
            print("knot exterior found, the knot is ")
            manifoldfound = True
    
    if manifoldfound == False:
        print("no knot exterior found")
        return False
    else:
        print("%d knot exterior found"%(len(exterior_list)))
        print(exterior_list)
        return True
    
print(check_knot_exterior())