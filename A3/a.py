import random

c = random.randint(1000,5000)
r = random.randint(1000,5000)

f = open('vect.txt','w')
for i in range(r):
    for j in range(c):
        n = random.uniform(-100, 100)
        f.write(str(n) + ' ')
    f.write('\n')

f.close()

