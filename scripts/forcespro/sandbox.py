
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from matplotlib import cm
import yaml
from matplotlib.animation import FuncAnimation
import os
import imageio


imagenames = []
Tf = 1.5
N = 40
Ts = Tf/N

for idx in range(355):
    imagenames.append('tempdump/f_'+str(idx)+'.jpg')


with imageio.get_writer("small_final.gif", mode='I', duration = Ts) as writer:
    for filename in imagenames:
        image = imageio.imread(filename)
        writer.append_data(image)


-3.2,
-3,
-2.4,
-1.8,
0,
0.6,
1.2,
1.4,
1.2,
0.6,
0,
-1.8,
-2.4,
-3,

-1.2,
-1.8,
-2.4,
-2.6,
2.6,
-2.4,
-1.8,
-1.2,
-0.6,
0,
0.2,
0.2,
0,
0.6
