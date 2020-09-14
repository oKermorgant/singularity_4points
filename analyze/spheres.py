#!/usr/bin/env python
import os
import sys
import subprocess
import yaml
from copy import deepcopy
import time

estims = ('dem','p4p','upnp','epnp')

def is_estim_folder(d):
    for e in estims:
        if d.startswith(e):
            return True
    return False

def pjoin(*f):
    return os.path.join(*f)

r_max = 1
result_dir=''

for arg in sys.argv[1:]:
    if arg.replace('.','').isdigit():
        r_max = float(arg)
    else:
        result_dir = arg
        
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
    
    

def run(cmd, wait = False):
    if wait:
        return subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    else:
        return subprocess.Popen(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    
# list simulated scenes
scenes = [pjoin(result_dir, d) for d in os.listdir(result_dir) if d.startswith('scene') and os.path.isdir(pjoin(result_dir, d))]

for scene in scenes:
    # find processed estimation methods
    methods = [d for d in os.listdir(scene) if is_estim_folder(d)]
    
    # find all performed simulations
    sims = set([d for m in methods for d in os.listdir(pjoin(scene, m)) if d.startswith('sphere')])
    
    for sim in sims:
        dst = os.path.join(scene, sim)
        if not os.path.exists(dst):
            os.mkdir(dst)
        
        # these folders have the same contents -> regroup them in a single file
        srcs = [m for m in methods if sim in os.listdir(pjoin(scene, m))]
        srcs.sort()
        
        # keep same ordering
        idx = 0
        last = srcs[0].split('_')[0]
        colors = []
        for k,src in enumerate(srcs):
            c = ''
            if src.endswith('vvs'):
                c = 's'
            elif src.endswith('lm'):
                c = 'o'
            
            if not src.startswith(last):
                last = src.split('_')[0]
                idx += 1
                
            colors.append('C{}-{}'.format(idx, c))
        
        srcs = [pjoin(scene,s,sim) for s in srcs]
        
        src_files = [f for f in os.listdir(srcs[0]) if f.endswith('.yaml') and not f.startswith('sphere')]        

        for src_file in src_files:
            
            dst_file = pjoin(dst, src_file)
            
            out = {'dataType': 'XY', 'xlabel': 'sphere radius [m]'}
            out['legend'] = []
            out['args'] = '--legendLoc out'
            out['lineType'] = colors
            
            init = False
            
            # concatenate all available simular simulations            
            for k,src in enumerate(srcs):
                
                with open(pjoin(src,src_file)) as f:
                    data = yaml.safe_load(f)
                    out['legend'].append(data['legend'][0])
                    if not init:
                        init = True
                        out['ylabel'] = data['ylabel']
                        out['data'] = [[float(val) for val in line] for line in data['data'] if line[0] <= r_max]
                    else:
                        for i,line in enumerate(data['data']):
                            if line[0] <= r_max:
                                out['data'][i] += line
            
            with open(dst_file, 'w') as f:
                yaml.safe_dump(out, f)
            print('-> ' + dst_file)
            run('log2plot {} --nodisplay'.format(dst_file))
        
        
