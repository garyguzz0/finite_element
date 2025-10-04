import numpy as np
import random
from matplotlib import pyplot as plt



class FiniteElement:
    def __init__(self,start,end,n):
        self.start = start
        self.end = end
        self.n = n
    
    def create_random_partition(self,start,end,n):
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
    
    def create_fixed_partition(self):
        start = 0
        end = 1
        n = 5

        bdrs = np.linspace(start,end,n+1)
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
        return np.linspace(interval[0],interval[-1],m+1,endpoint=True)
    
    def derivative(self,phi,omega):
        dy = np.gradient(phi,omega)
        return dy[:-1]

    def create_shape_functions(self,partition):
        Phi = []
        Omega = []
        Dx = []
        
        sf1 = self.shape_function_1
        sf2 = self.shape_function_2
        for i in range(self.n):
            interval = partition[i]
            element_domain = self.create_element_domain(interval)
            dx = element_domain[1]-element_domain[0]
            phi = []
            
            phi.append([sf1(x,interval) for x in element_domain])
            phi.append([sf2(x,interval) for x in element_domain])
            
            Phi.append(np.array(phi))
            Omega.append(element_domain)
            Dx.append(dx)
            
        return [Phi,Omega,Dx]
    
    def plot_shape_functions(self,Phi,Omega):
        colors = plt.cm.Set2(np.linspace(0,1,self.n))
        for i in range(len(Phi[0])):
            for j in range(self.n):
                plt.plot(Omega[j],Phi[j][i,:],color=colors[j])
        plt.show()
    
    def compute_stiffness_entry(self,phi_1,phi_2,omega,dx):
        dphi_1 = self.derivative(phi_1,omega)
        dphi_2 = self.derivative(phi_2,omega)
        return np.dot(dphi_1,dphi_2)*dx
    
    def compute_stiffness_matrix(self,phi,omega,dx):
        k = np.zeros((2,2))
        for i in range(k.shape[0]):
            for j in range(k.shape[1]):
                k[i,j] = self.compute_stiffness_entry(phi[i,:],phi[j,:],omega,dx)
        return k
    
    def compute_stiffness_matrices(self,Phi,Omega,Dx):
        K_e = []
        for i in range(self.n):
            K_e.append(self.compute_stiffness_matrix(Phi[i],Omega[i],Dx[i]))
        return K_e
    
    def assemble_global_stiffness_matrix(self,K_e):
        K = np.zeros((self.n+1,self.n+1))
        for i in range(self.n):
            K[i:i+2, i:i+2] += K_e[i]
        return K
    
def main():
    start = 0
    end = 1
    n=5

    fem = FiniteElement(start,end,n)
    partition = fem.create_random_partition(start, end, n)
    # partition = fem.create_fixed_partition()

    [Phi,Omega,Dx] = fem.create_shape_functions(partition)
    
    fem.plot_shape_functions(Phi,Omega)

    K_e = fem.compute_stiffness_matrices(Phi,Omega,Dx)
    K = fem.assemble_global_stiffness_matrix(K_e)
    
    return K

if __name__ == "__main__":
    K = main()
    
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

