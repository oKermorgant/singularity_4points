#!/usr/bin/env python
import os, shutil
import sys
import subprocess
import yaml
from copy import deepcopy
import time

estims = [('dem','dem_vvs', 'dem_lm'),('p4p',),('upnp','upnp_vvs','upnp_lm')]
exp_dir = 'rose_noNoise'
trans_error = 'te.yaml'
pose_gt = 'cMs.yaml'
pose_raw = 'iMs.yaml'
pose_final = 'ceMs.yaml'


def is_estim_folder(d):
    for e in estims:
        if d.startswith(e):
            return True
    return False

def pjoin(*f):
    return os.path.join(*f)

result_dir=''
        
if result_dir == '':
    with open('../config.yaml') as f:
        config = yaml.safe_load(f)
    result_dir = config['dataPath']
    
if not os.path.exists(result_dir):
    # run from the code folder
    result_dir = '../results'
    
if not os.path.exists(result_dir):
    print('Folder {} does not exist'.format(result_dir))
    sys.exit(0)
    
result_dir += 'scene1/'
    
def run(cmd, wait = False):
    if wait:
        return subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    else:
        return subprocess.Popen(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    
# group all similar errors (to copy to correct dir)
for estim in estims:    
    files = [pjoin(result_dir, e, 'rose_noNoise/te.yaml') for e in estim]
    run('log2plot {} --nodisplay -g'.format(' '.join(files)))
    
    
for estim in estims:
   srcs = [pjoin(result_dir, e, 'rose_noNoise/ceMs.yaml') for e in estim]
   srcs.insert(0, pjoin(result_dir, estim[0], 'rose_noNoise/cMs.yaml'))
   
   dst = []
   for i,src in enumerate(srcs):
        with open(src) as f:
            data = yaml.safe_load(f)
        data['lineType'] = ['C{}'.format(i==0 and 3 or i-1)]
        dst.append(pjoin(result_dir, 'rose_noNoise/{}_{}.yaml'.format(estim[0],i)))
        with open(dst[-1], 'w') as f:
            yaml.safe_dump(data, f)
       
   run('log2plot {}'.format(' '.join(dst)))
    
    
    
       
    
