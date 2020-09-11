#!/usr/bin/env python
import os
import sys
import subprocess
import yaml
from copy import deepcopy
import time

estims = ('dem','p4p','upnp')
ref = ('','_VVS', '_LM')
r_max = len(sys.argv) == 2 and float(sys.argv[1]) or 1

def run(cmd, wait = False):
    if wait:
        return subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
    else:
        return subprocess.Popen(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
for config in range(1,11):
    basedir = 'config{}'.format(config)

    now = time.time()
    
    for noise in (0, 0.01, 0.001, 0.0001, 0.00001):
        
        noise_s = noise and f'{noise:.7f}'.strip('0')[1:] or ''
        #os.chdir(os.path.abspath(os.path.dirname(__file__)))
        
        if noise_s:
            dest_dir = '{}/sphere_{}'.format(basedir, noise_s)
        else:
            dest_dir = '{}/sphere_perfect'.format(basedir)
            
            
        if not os.path.exists(dest_dir):
            os.mkdir(dest_dir)
        
        dest_concat = dest_dir + '/summary.pdf'

        def src_dir(e):
            if noise:
                return '{}/{}/noise_{}/'.format(basedir, e, noise_s)
            return '{}/{}/'.format(basedir,e)
        
        if not os.path.exists(src_dir(estims[0])):
            print(dirname(estims[0]) + ' does not exists')
            continue

        descriptions = [fi for fi in os.listdir(src_dir(estims[0])) if fi.startswith('sphere') and fi.endswith('_err.yaml') and 'VVS' not in fi and 'LM' not in fi]
        descriptions = [des.replace('sphere_', '').replace('_err.yaml', '') for des in descriptions]

        dests = []

        for des in descriptions:
            
            out = {'dataType': 'XY', 'xlabel': 'sphere radius [m]'}
            out['legend'] = []
            
            init = False

            for e in estims:
                for r in ref:
                    src_file = src_dir(e) + 'sphere{}_{}_err.yaml'.format(r, des)
                    with open(src_file) as f:
                        data = yaml.safe_load(f)
                    #print(src_file)
                    #print('  -> dim = {}'.format(len(data['data'])))

                    out['legend'].append(data['legend'][0].replace('P3P', 'P4P'))
                    
                    if not init:
                        init = True
                        out['ylabel'] = data['ylabel']
                        out['data'] = [[float(val) for val in line] for line in data['data'] if line[0] <= r_max]
                        m = max(val[1] for val in data['data'])
                    else:
                        m = max(m, max(val[1] for val in data['data']))
                        for i,line in enumerate(data['data']):
                            if line[0] <= r_max:
                                out['data'][i] += line
                    
                    # also create 3D sphere plot
                    #if des == descriptions[0]:
                    #    src_file = dirname(e) + 'sphere{}_3D.yaml'.format(r)
                    #    run('log2plot {} --nodisplay'.format(src_file))
                            
            out['args'] = '--legendLoc out'
                
            dest = '{}/{}.yaml'.format(dest_dir, des)
            pdf = dest.replace('.yaml', '.pdf')
            if os.path.exists(pdf):
                os.remove(pdf)
            dests.append(pdf)
            with open(dest, 'w') as f:
                yaml.safe_dump(out, f)
            #os.system('rosrun log2plot plot {} --nodisplay &'.format(dest))
            run('log2plot {} --nodisplay'.format(dest))
            
        # join all pdf's
        print('   Creating from ' + ' '.join(dests))

        for d in dests:
            if not os.path.exists(d):
                time.sleep(1)
        time.sleep(1)
        out = run('pdftk {} cat output {}'.format(' '.join(dests), dest_concat), True)
    
