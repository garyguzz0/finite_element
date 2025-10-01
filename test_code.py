import numpy as np

# Material and geometric properties
E = 200e9  # Young's modulus (Pa)
A = 0.001  # Cross-sectional area (m^2)
L = 1.0    # Total length of the bar (m)
num_elements = 2  # Number of elements

# Applied force
P = 1000   # Axial force (N)

# Calculate element length
element_length = L / num_elements

# Initialize global stiffness matrix and force vector
num_nodes = num_elements + 1
K_global = np.zeros((num_nodes, num_nodes))
F_global = np.zeros(num_nodes)

# Assemble global stiffness matrix
for i in range(num_elements):
    k_element = (A * E / element_length) * np.array([[1, -1], [-1, 1]])
    
    # Map local DOFs to global DOFs
    global_dofs = [i, i + 1]
    
    for row_idx, global_row in enumerate(global_dofs):
        for col_idx, global_col in enumerate(global_dofs):
            K_global[global_row, global_col] += k_element[row_idx, col_idx]

# Apply boundary conditions (fixed at node 0)
# Remove row and column corresponding to fixed DOF
fixed_dof = 0
K_reduced = np.delete(np.delete(K_global, fixed_dof, axis=0), fixed_dof, axis=1)

# Apply applied load (at node 'num_nodes - 1')
F_global[num_nodes - 1] = P
F_reduced = np.delete(F_global, fixed_dof)

# Solve for displacements
U_reduced = np.linalg.solve(K_reduced, F_reduced)

# Reconstruct full displacement vector
U_global = np.insert(U_reduced, fixed_dof, 0)

print("Global Stiffness Matrix:")
print(K_global)
print("\nGlobal Force Vector:")
print(F_global)
print("\nDisplacement at each node:")
print(U_global)

# Calculate strain and stress in each element (optional)
for i in range(num_elements):
    u1 = U_global[i]
    u2 = U_global[i+1]
    
    strain = (u2 - u1) / element_length
    stress = E * strain
    print(f"\nElement {i+1}:")
    print(f"Strain: {strain:.6e}")
    print(f"Stress: {stress:.2f} Pa")
