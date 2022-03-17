indptr = []
index = []
levelo = []

with open('dummy/level_offset.txt', 'r') as f:
    for line in f:
        levelo.append(int(line.strip('\n')))
with open('dummy/indptr.txt', 'r') as f:
    for line in f:
        indptr.append(int(line.strip('\n')))
with open('dummy/index.txt', 'r') as f:
    for line in f:
        index.append(int(line.strip('\n')))

print('indptr')
print(indptr)

print('index')
print(index)

print('levelo')
print(levelo)