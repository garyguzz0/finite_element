import numpy as np
import random
from matplotlib import pyplot as plt
import scipy


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
        
    # N_0(x)
    def shape_function_1(self,x,interval):
        start = interval[0]
        end = interval[-1]
        
        width = end-start
        m = 1/width
        return -m*(x-end)
         
    # N_1(x)
    def shape_function_2(self,x,interval):
        start = interval[0]
        end = interval[-1]
        
        width = end-start
        m = 1/width
        return m*(x-start)
    
    def create_element_domain(self,interval,m=100):
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
    
    def f(self,x):
        return 1
    
    def create_f_vector(self,omega):
        f = self.f
        return np.array([f(x) for x in omega])
        
    def compute_element_vector(self,f_vector,phi,dx):
        return np.dot(-f_vector,phi)*dx
    
    def compute_element_vectors(self,Phi,Omega,Dx):
        F_e = []
        create_f_vector = self.create_f_vector
        compute_element_vector = self.compute_element_vector
        
        for i in range(self.n):
            F_e.append([])
            for j in range(2):
                f_vector = create_f_vector(Omega[i])
                element_vector = compute_element_vector(f_vector,Phi[i][j,:],Dx[i])
                
                F_e[i].append(element_vector)
        return F_e
    
    def assemble_global_element_vector(self,F_e):
        F = np.zeros(self.n+1)
        for i in range(self.n):
            F[i:i+2] += F_e[i]
        return F
    
    # cN(xi)
    def curly_N(self,xi,node):
        conj = -2*node+1
        return 1/2-conj*1/2*xi
    
    # phi_e(xi) = x
    def phi_e(self,xi,a,b):
        cN = self.curly_N
        x = cN(xi,0)*a + cN(xi,1)*b
        return x
    
    # phi_e : omega_ref -> omega
    def x_vector(self,omega_ref,a,b):
        phi_e = self.phi_e
        omega = np.array([phi_e(xi,a,b) for xi in omega_ref])
        return omega
    
    # phi_e_inverse(x) = xi
    def phi_e_inverse(self,x,a,b):
        return (1/2*(a+b)-x)/(1/2*(a-b))
    
    # phi_e_inverse : omega -> omega_ref
    def xi_vector(self,omega):
        phi_e_inverse = self.phi_e_inverse
        a = omega[0]
        b = omega[-1]
        omega_ref = np.array([phi_e_inverse(x,a,b) for x in omega])
        return omega_ref
    
    def compute_partial_N_partial_x(self,omega_ref,a,b,node):
        curly_N = self.curly_N
        phi_e = self.phi_e
        
        partial_cN_partial_xi = scipy.differentiate.derivative(curly_N, omega_ref,args=node).df
        partial_x_partial_xi = scipy.differentiate.derivative(phi_e, omega_ref, args=(a,b)).df
        
        return partial_cN_partial_xi * (1/partial_x_partial_xi)
    
    def compute_K_ij_e(self,omega_ref,a,b,i,j):
        dxi = np.gradient(omega_ref)
        partial_cNi_parital_xi = scipy.differentiate.derivative(self.curly_N,omega_ref,args=i).df
        partial_cNj_partial_xi = scipy.differentiate.derivative(self.curly_N,omega_ref,args=j).df
        
        partial_Ni_partial_x = self.compute_partial_N_partial_x(omega_ref,a,b,i)
        partial_Nj_partial_x = self.compute_partial_N_partial_x(omega_ref,a,b,j)
        partial_x_partial_xi = scipy.differentiate.derivative(self.phi_e, omega_ref, args=(a,b)).df
        
        integrand = (partial_Ni_partial_x)*(partial_Nj_partial_x)*(partial_x_partial_xi)*dxi
        return scipy.integrate.trapezoid(integrand)
    
    def compute_F_i_e(self,omega_ref,a,b,i):
        dxi = np.gradient(omega_ref)
        partial_x_partial_xi = scipy.differentiate.derivative(self.phi_e, omega_ref, args=(a,b)).df
        curly_N = self.curly_N
        curly_N_vector = np.array([curly_N(xi,i) for xi in omega_ref])
        x_vector = self.x_vector(omega_ref,a,b)
        f = self.f
        f_vector = np.array([f(x) for x in x_vector])
        
        integrand = -f_vector*curly_N_vector*partial_x_partial_xi*dxi
        return scipy.integrate.trapezoid(integrand)
    
    def get_gauss_points_and_weights(self,order):
        if order <= 1:
            return {"points":0,"weights":2}
        elif 1 < order <= 3:
            return {"points":[-1/np.power(3,1/2),1/np.power(3,1/2)],"weights":[1,1]}
        elif 3 < order <=5:
            return {"points":[-np.power(3/5,1/2), 0, np.power(3/5,1/2)],"weights":[5/9, 8/9, 5/9]}
        
    def evaluate_integrand_K(self,xi,i,j,a,b):
        cN_i = self.curly_N
        cN_j = self.curly_N
        phi_e = self.phi_e
        
        d_cN_i = scipy.differentiate.derivative(cN_i,xi,args=i).df
        d_cN_j = scipy.differentiate.derivative(cN_j,xi,args=j).df
        
        d_phi_e = scipy.differentiate.derivative(phi_e,xi,args=(a,b)).df
        
        return d_cN_i*(1/d_phi_e)*d_cN_j*(1/d_phi_e)*d_phi_e      
    
    def gaussian_quadrature(self,integrand,order):
        gpw = self.get_gauss_points_and_weights(order)
        
        return np.dot(integrand,gpw['weights'])
    
    
        
    # def transform_reference_domain(self,x,)
    
    # def compute_reference_element(self):
        
        
    
    
def main():
    start = 0
    end = 1
    n=5

    fem = FiniteElement(start,end,n)
    # partition = fem.create_random_partition(start, end, n)
    partition = fem.create_fixed_partition()

    [Phi,Omega,Dx] = fem.create_shape_functions(partition)
    
    omega = Omega[0]

    omega_ref = fem.xi_vector(omega)

    partial_N_partial_x = fem.compute_partial_N_partial_x(omega_ref,omega[0],omega[-1],0)
    
    K_ij_e = fem.compute_K_ij_e(omega_ref,omega[0],omega[-1],0,1)
    F_i_e = fem.compute_F_i_e(omega_ref, omega[0], omega[-1],0)

    gpw = fem.get_gauss_points_and_weights(1) 
    
    integrand = fem.evaluate_integrand_K(0,0,1,0,.2)
    
    gq = fem.gaussian_quadrature(integrand,1)
    
    # fem.plot_shape_functions(Phi,Omega)

    K_e = fem.compute_stiffness_matrices(Phi,Omega,Dx)
    K = fem.assemble_global_stiffness_matrix(K_e)
    
    F_e = fem.compute_element_vectors(Phi,Omega,Dx)
    F = fem.assemble_global_element_vector(F_e)
    
    return K, F

if __name__ == "__main__":
    [K,F] = main()
    # F = np.array([.1,.2,.2,.2,.2,.1]) This is what F gives when computed via integrals (as opposed to Riemann sums), but even using this in linalg.solve, it's still wrong
    # U = np.linalg.solve(K,F) #This is where it goes wrong
    # U = scipy.linalg.solve(K,-F)
    
    # plt.plot(np.linspace(0,1,6),U)
    # plt.show()

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

