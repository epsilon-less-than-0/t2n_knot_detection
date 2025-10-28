
from snappy import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
b=BraidGroup(7)([1, 2, 1, -2, 3, 2, 1, -2, -3, 4, 3, 2, 1, -2, -3, -4, 5, 4, 3, 2, -3, -4, -5, 6, 5, 4, 3, -4, -5, -6])
l = list(b.Tietze())
M = Manifold('Braid%s(2,0)(0,0)'%(l))

homologyZZ = lambda Y: Y.homology().elementary_divisors() == [0]
coverlist = list(filter(homologyZZ, M.covers(2)))

if len(coverlist) == 1:
    print("single manifold in coverlist so double cover found")

G = coverlist[0]

L = G.filled_triangulation().exterior_to_link(verbose=True)

print(L.knot_floer_homology())
 