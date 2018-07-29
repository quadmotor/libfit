
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def tokenizer(fname):
    with open(fname) as f:
        chunk = []
        for line in f:
            if 'HEAD'in line:
                continue
            if 'END' in line:
                yield chunk
                chunk = []
                continue
            chunk.append(line)


xs = [np.loadtxt(A) for A in tokenizer('output/x_cp.dat')]
ys = [np.loadtxt(A) for A in tokenizer('output/y_cp.dat')]
zs = [np.loadtxt(A) for A in tokenizer('output/z_cp.dat')]
for i in range(len(xs)):
       # ax.plot_wireframe(xs[i], ys[i], zs[i], rstride=1, cstride=1, color='g')
	ax.scatter(xs[i], ys[i], zs[i], color='g')


xs = [np.loadtxt(A) for A in tokenizer('output/x_in.dat')]
ys = [np.loadtxt(A) for A in tokenizer('output/y_in.dat')]
zs = [np.loadtxt(A) for A in tokenizer('output/z_in.dat')]
#for i in range(len(xs)):
#	ax.plot_surface(xs[i], ys[i], zs[i], rstride=1, cstride=1, color='r')
	
	
xs = [np.loadtxt(A) for A in tokenizer('output/x_out.dat')]
ys = [np.loadtxt(A) for A in tokenizer('output/y_out.dat')]
zs = [np.loadtxt(A) for A in tokenizer('output/z_out.dat')]
for i in range(len(xs)):
	ax.plot_wireframe(xs[i], ys[i], zs[i], rstride=1, cstride=1, color='b')
	#ax.scatter(xs[i], ys[i], zs[i], color='b')

ax.set_xlabel('x (u, v)')
ax.set_ylabel('y (u, v)')
ax.set_zlabel('z (u, v)')
plt.show()
