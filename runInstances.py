import os
import re
import sys

pasta_instancias = 'instances'
comando_exec = 'python3 exact\\(stages\\)/main.py'
pasta_resultados = 'results'
maxtam = int(sys.argv[-1]) if len(sys.argv) >= 2 else 500

os.system('ulimit -Sv 4000000') # Limite máximo de RAM = 4GB
os.makedirs(pasta_resultados, exist_ok=True)

for f in os.listdir(pasta_instancias):
    if not f.endswith('.tsp') or int(re.findall('\\d+', f)[0]) >= maxtam:
        continue
    outf = os.path.join(pasta_resultados, f)
    if os.path.isfile(outf):
        continue
    s = '{} {} >> {}'.format(comando_exec, os.path.join(pasta_instancias, f), outf)
    #print(s)
    os.system(s)