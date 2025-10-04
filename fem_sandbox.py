import numpy as np
import random
from matplotlib import pyplot as plt



class FiniteElement:
    def __init__(self,start,end,n):
        self.start = start
        self.end = end
        self.n = n
    
    def create_partition(self,start,end,n):
        start = self.start
        end = self.end
        n = self.n
        
        bdrs = set()
        while len(bdrs) < n - 1:
            bdrs.add(random.uniform(start,end))
        bdrs = sorted(list(bdrs))
        bdrs = [start] + bdrs + [end]
        
        partition = []
        for i in range(len(bdrs)-1):
            partition.append((bdrs[i],bdrs[i+1]))
        return partition
        
    def visualize_interval_partitions(self, partitions, original_start, original_end, title="Interval Partitions"):
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
        
    def shape_function_1(self,x,interval):
        start = interval[0]
        end = interval[-1]
        
        width = end-start
        m = 1/width
       
        return -m*(x-end)
            
    def shape_function_2(self,x,interval):
        start = interval[0]
        end = interval[-1]
        
        width = end-start
        m = 1/width
        
        return m*(x-start)
    
    def create_element_domain(self,interval,m=50):
        return np.linspace(interval[0],interval[-1],m,endpoint=True)

    def create_shape_functions(self,partition):
        Phi = []
        
        sf1 = self.shape_function_1
        sf2 = self.shape_function_2
        for i in range(self.n):
            interval = partition[i]
            element_domain = self.create_element_domain(interval)
            phi = []
            
            phi[0] = [sf1(x,interval) for x in element_domain]
            phi[1] = [sf2(x,interval) for x in element_domain]
            
            Phi.append(phi)
        
        return Phi
            

start = 0
end = 1
n = 8

fem = FiniteElement(start,end,n)
partition = fem.create_partition(start, end, n)

Phi = fem.create_shape_functions(partition)

