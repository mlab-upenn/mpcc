
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from matplotlib import cm
import yaml
from matplotlib.animation import FuncAnimation
import os
import imageio

solverparams = "parameters/solverparams.yaml"
#load global constant model parameters
with open(solverparams) as file:
    params = yaml.load(file, Loader= yaml.FullLoader)

print(params["Qc"])
