import numpy as np
import random
from matplotlib import pyplot as plt



class FiniteElement:
    def __init__(self,start,end,n):



def create_partition(start,end,n):
    bdrs = set()
    while len(bdrs) < n - 1:
        bdrs.add(random.uniform(start,end))
    bdrs = sorted(list(bdrs))
    bdrs = [start] + bdrs + [end]
    
    partition = []
    for i in range(len(bdrs)-1):
        partition.append((bdrs[i],bdrs[i+1]))
    return partition
    
def visualize_interval_partitions(partitions, original_start, original_end, title="Interval Partitions"):
    fig, ax = plt.subplots(figsize=(10, 2))
    
    ax.hlines(0, original_start, original_end, color='gray', linewidth=1, linestyle='--', zorder=0)
    colors = plt.cm.viridis(np.linspace(0, 1, len(partitions))) 
    
    for i, (start, end) in enumerate(partitions):
        ax.hlines(0, start, end, color=colors[i], linewidth=5, zorder=1) 
        ax.plot([start, end], [0, 0], marker='|', markersize=10, color=colors[i], zorder=2)

    ax.set_ylim(-0.1, 0.2)
    ax.set_xlim(original_start, original_end)
    ax.yaxis.set_visible(False)
    ax.tick_params(axis='x', length=6, width=1.5, color='black')

    plt.title(title)
    plt.show()
    
def shape_function(start,end):
    width = end-start
    phi_0 = 
    
def create_shape_functions(partition):
    n = len(partition)
    
    for i = 
    



start = 0
end = 1
n = 8

partition = create_partition(start,end,n)
visualize_interval_partitions(partition,start,end)