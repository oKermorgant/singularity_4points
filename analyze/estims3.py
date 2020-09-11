from pylab import *
import os
import yaml
from copy import deepcopy
from matplotlib import rc, animation, rcParams
from matplotlib.font_manager import FontProperties

# LateX
rc('font', family='sans-serif')
rc('text', usetex = True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rcParams['text.latex.preamble']=[
    r"\usepackage{amsmath}",
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath']

def fuse(l1, l2):
    ret = [0] * (len(l1) + len(l2))
    ret[::2] = l1
    ret[1::2] = l2
    return ret

leg = (('dem', 'DeMenthon'),('p4p', 'P4P'), ('epnp', 'EPnP'))

class Error:
    def __init__(self):
        self.vals = []
        
    def read(self, filename): 
        
        with open(filename) as f:
            error = yaml.load(f)['data']
        if len(error[0]) == 6:
            # pose error -> keep translation norm
            error = [norm(line[:3]) for line in error]

        self.vals.append(error)
        
    def get(self, des):
        meth = {'max': max, 'min': min, 'mean': mean, 'median': median, 'std': std}
        return [meth[des](e) for e in self.vals]

class Estimator:
    def __init__(self, method):
        self.method = method
        for m,l in leg:
            if m == method:
                self.legend = l
                break
        
        # load files
        self.t = Error()
        self.t_vvs = Error()
        self.z = Error()
        self.z_vvs = Error()
        
    def read(self):
        
        for c in range(10):
            
            path = 'config{}/{}/rose_'.format(c+1, self.method)            
            self.t.read(path + 'cMce.yaml')
            self.z.read(path + 'Z.yaml')

            path = 'config{}/{}/rose_VVS_'.format(c+1, self.method)            
            self.t_vvs.read(path + 'cMce.yaml')
            self.z_vvs.read(path + 'Z.yaml')
            
estim = [Estimator(m[0]) for m in leg]

#if not os.path.exists('estim_summary.yaml'):
if not False:
    summary = {}
    
    for e in estim:
        print('Reading ' + e.method)
        e.read()
        #scores = {}

        #for field in ('t','z','t_vvs','z_vvs'):
            #err = getattr(e, field)
            #scores[field] = {'vals': err.vals, 'min': err.vmin}                
        #summary[e.method] = deepcopy(scores)
    #with open('estim_summary.yaml', 'w') as f:
        #yaml.dump(summary, f, default_flow_style=False)        
#else:
    #with open('estim_summary.yaml') as f:
        #summary = yaml.load(f)
        #for e in estim:
            #scores = deepcopy(summary[e.method])
            #for field in ('t','z','t_vvs','z_vvs'):
                #getattr(e, field).vals = scores[field]['vals']
                #getattr(e, field).vmin = scores[field]['min']

close('all')

N = len(estim)
ind = arange(10.)
width = .8/(2*N)
#c = ([0,.6,0],[0,1,0],[.6,0,.5],[0,0,1],[.8,.6,0],[1,0,0])
size = 15

for out in ('max', 'mean','median'):
    use_max = out == 'max'
    
    ft = figure(figsize=(8,4))
    pt = ft.gca()
    pt.set_axisbelow(True)


    for i,e in enumerate(estim):
        offset = (2*i-N/2-1.5)*width
        k = 0
        for t,lab in ((e.t,''), (e.t_vvs, ' +VVS')):
            pt.bar(ind+offset, t.get(out), width, label=e.legend+lab)
            
            if out == 'mean':
                m = t.get('min')
                M = t.get('max')
                for j,p in enumerate(ind+offset):
                    pt.plot([p,p],[m[j],M[j]], 'k', linewidth=2)               
            
            k += 1            
            offset += width

    pt.legend(ncol=N, prop= FontProperties(size=size-2))

    pt.set_ylabel('{} translation error [m]'.format(out.title()), size = size+2)

    for p in [pt]:
        p.set_xticks(ind)
        p.set_xticklabels(['Ex. {}'.format(int(i+1)) for i in ind])
        
        # double y ticks
        dy = p.get_yticks()[1] - p.get_yticks()[0]
        ym,yM = p.get_ylim()
        yM = ym+ceil((yM-ym)*1.15/dy)*dy
        p.set_ylim(ym,yM)
        p.set_yticks(arange(ym, yM+dy,dy))
        #m = p.get_yticks()[1]
        #M = m*ceil(p.axis()[-1]/m)
        
        p.grid(which='both', axis='y')
                    
        for ti in p.get_xticklabels():
            ti.set_fontsize(size)
        for ti in p.get_yticklabels():
            ti.set_fontsize(size)

    ft.tight_layout()
    ft.savefig('err_t_{}.pdf'.format(out))
    
    draw()
show()




'''
fz = figure()
pz = fz.gca()

for c in range(10):
    errors = fuse([e.t.vals[c] for e in estim],[e.t_vvs.vals[c] for e in estim])
    pt.bar(ind, errors, width, bottom=bottomt)
    for i in range(2*N):
        bottomt[i] += errors[i]
    errors = fuse([e.z.vals[c] for e in estim],[e.z_vvs.vals[c] for e in estim])
    pz.bar(ind, errors, width, bottom=bottomz)
    for i in range(2*N):
        bottomz[i] += errors[i]
        
size = 15

pt.set_ylabel('Max translation estimation error [m]', size = size+2)
pz.set_ylabel('Max depth estimation error [\%]', size = size+2)

for p in [pt, pz]:
    p.set_xticks(ind)
    p.set_xticklabels(fuse([leg[e.method] for e in estim], ['+VVS' for _ in estim]))
    
    # double y ticks
    m = p.get_yticks()[1]
    M = m*ceil(p.axis()[-1]/m)
    #p.set_yticks(arange(0, M+m/2, m/2))
    
    p.grid(which='both', axis='y')
                 
    for ti in p.get_xticklabels():
        ti.set_fontsize(size)
    for ti in p.get_yticklabels():
        ti.set_fontsize(size)

ft.tight_layout()
fz.tight_layout()
'''

draw()
show()

ft.savefig('err_t.pdf')

#fz.savefig('err_z.pdf')
