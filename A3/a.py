import random

c = random.randint(10,100)
r = random.randint(10,100)

f = open('vect.txt','w')
for i in range(r):
    for j in range(c):
        n = random.uniform(-100, 100)
        f.write(str(n) + ' ')
    f.write('\n')

f.close()