from pylab import *
import yaml
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.font_manager import FontProperties
import os

legendSize = FontProperties(size=18)
markersize = 10
labelsize = 18



with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
    
desired = data['fixedObject1']['nodes'][0]

Z = array([[v == 'nan' and nan or float(v) for v in line] for line in data['data']])

# median filter step for remaining outliers
Zm = Z[1:-2,1:-2].copy()
m = amax(Zm)
for i in range(1, Z.shape[0]-2):
    for j in range(1, Z.shape[0]-2):
        Zm[i-1,j-1] = median(Z[i-1:i+2,j-1:j+2])/m
        Zm[i-1,j-1] = -1./(.15+(Zm[i-1,j-1]))
side = data['args']
print('Steps: ',Z.shape[0])

close('all')
fig = plt.figure(figsize=(8.3,8))
ax = fig.gca()
axis('equal')

alpha = .6

pos = imshow(Zm.T, extent=(-side,side,-side,side), origin='lower')
#,vmax=alpha*mean(Z) + (1-alpha)*amax(Z))
#fig.colorbar(pos, ax=ax)

if 'o/no_estim' not in os.path.abspath(sys.argv[1]):
    ax.plot([0], [0], 'rD', label='Singular point', markersize=markersize+2)
ax.plot([desired[0]], [desired[1]], 'ys', label='Desired position', markersize=markersize)

axis((-side,side,-side,side))
xlabel('$\Delta X$ [m]', size = labelsize)
ylabel('$\Delta Y$ [m]', size = labelsize)

legend(prop=legendSize)

for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontsize(labelsize)
    
tight_layout()

show()

savefig(sys.argv[1].replace('.yaml','.pdf'))
