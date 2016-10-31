import time as T
import random as R

print "Running method 1"

R.seed(32)
hopHistory = []
t0 = T.time()
for i in range(10000):
    randNo = R.randint(0, 50000)
    hopHistory.append(randNo)
    hopHistory = list(set(hopHistory))
t1 = T.time()


print "Running method 2"

R.seed(32)
hopHistory2 = []
t2 = T.time()
for i in range(10000):
    randNo = R.randint(0, 50000)
    if randNo not in hopHistory2:
        hopHistory2.append(randNo)
t3 = T.time()


print "Running method 3"

R.seed(32)
hopHistory3 = []
t4 = T.time()
for i in range(10000):
    randNo = R.randint(0, 50000)
    hopHistory3.append(randNo)
hopHistory3 = list(set(hopHistory3))
t5 = T.time()


print "RESULTS"
# if (hopHistory != hopHistory2) or (hopHistory != hopHistory3) or (hopHistory2 != hopHistory3):
#     raise SystemError('ISSUE')
print "Method 1 (set every hop) =", t1-t0
print "Method 2 (check every hop then append) =", t3-t2
print "Method 3 (set after carrier completion) =", t5-t4
