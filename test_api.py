import subprocess as sub
import time
import shutil
import os
from_file = "POTCAR"
to_file = "TTTT/POTCAR"
a = time.time()
for i in range(10000):
    with open(from_file) as file:
        dd = file.readline()
b = time.time()
print(b-a)

a = time.time()
for i in range(100):
    os.popen('sed -n {}p {}'.format(500, from_file)).read()
b = time.time()
print(b-a)