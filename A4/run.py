import subprocess
from sys import argv

folder = argv[1]
f = open(folder+"/arg.txt")
args = []
for i in range(3):
    l = float(f.readline().split('=')[1].strip())
    args.append(l)
subprocess.run(f'make test data={folder+"/data_image_a.txt"} query={folder+"/query_image_a.txt"} t1={args[1]} t2={args[2]} n={args[0]}',shell=True)