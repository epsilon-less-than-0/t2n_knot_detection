import time
start_time = time.time()
b = BraidGroup(7)([1, 2, 1, -2, 3, 2, 1, -2, -3, 4, 3, 2, 1, -2, -3, -4, 5, 4, 3, 2, -3, -4, -5, 6, 5, 4, 3, -4, -5, -6])

still_successful = False
power = 4
while still_successful == False: 
  print("checking %dth power"%(power))
  ap = (b**power).alexander_polynomial()
  if abs(ap.coefficients()[0]) == 1:
    print("*** %dth power of braid is fibered so might be exchangeable"%(power))
    power += 1
    print("--- %s seconds ---" % (time.time() - start_time))
  else: 
    still_successful = True
    print("%dth power of braid is not fibered"%(power))
    print("--- %s seconds ---" % (time.time() - start_time))  