import numpy as np
import random
from matplotlib import pyplot as plt


start = 0
end = 1
num_partitions = 8



def create_partition(start,end,num_partitions):
    bdrs = set()
    while len(bdrs) < num_partitions - 1:
        bdrs.add(random.uniform(start,end))
    bdrs = sorted(list(bdrs))
    bdrs = [start] + bdrs + [end]
    
    prtn = []
    for i in range(len(bdrs)-1):
        prtn.append((bdrs[i],bdrs[i+1]))
    return prtn
    
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
    

prtn = create_partition(start,end,num_partitions)
visualize_interval_partitions(prtn,start,end)