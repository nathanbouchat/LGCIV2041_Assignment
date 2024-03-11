# This code was originally prepared for the course LGCIV1023 (Stabilité des Constructions), on Decembre 3 2019
# Igor Bouckaert, João Pacheco de Almeida, Martin Steinmetz
# Adapted for the course "Nonlinear Response Analysis", Rose, Pavia, Italy, in May 2021
# Adapted for the course LGCIV2041 (Numerical Analysis of Civil Engineering Structures) by João Pacheco de Almeida on February 1 2023
# Application of the displacement method to solve frame structures

"""### 0: INTRODUCTION

This code written for Python / Google Colaboraty / Jupyter notebook aims to go step by step through a numerical solution of a frame structure with the stiffness or displacement method.

To do this, simply run each cell of code one after the other and check what each one does.

#### IMPORTANT

Before starting, it is also necessary to run all the cells signaled with #DEFINITIONS within this section header, which define the different functions and libraries that will be used in the rest of the code.
"""

# DEFINITIONS
# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

# DEFINITIONS
# Display of the element undeformed configuration based on the connectivity (incidence) matrix

def PlotUndeformed(coord, connect) : 

    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(Coord[0], Coord[1], 'ro')
    
    for i in range(len(coord[1])) : 
        plt.annotate(str(i), (coord[0][i], coord[1][i]))
    
    plt.show()

# DEFINITIONS
# Display of the deformed configuration of the element

def PlotDeformed(coord, connect, displ, scale) : 
    
    coord_def = coord + displ * scale
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
        plt.plot([coord_def[0][connect[i][0]] , coord_def[0][connect[i][1]]] , 
                [coord_def[1][connect[i][0]] , coord_def[1][connect[i][1]]] , 
                'r-', linewidth=0.5)
    
    plt.plot(coord[0], coord[1], 'ro')
    
    for i in range(len(coord[1])) : 
        plt.annotate(str(i), (coord[0][i], coord[1][i]))
    
    plt.show()


def PlotTransversalDisplacement(coord, connect, value, coord_exact, connect_exact, value_exact, n_elem):
    for i in range(len(connect_exact)) : 
        plt.plot( [coord_exact[0][connect_exact[i][0]] , coord_exact[0][connect_exact[i][1]]] , 
                [value_exact[connect_exact[i][0]] , value_exact[connect_exact[i][1]]] , 
                'b-', linewidth=0.5)
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value[connect[i][0]] , value[connect[i][1]]] , 
                'r--', linewidth=0.5)
    plt.plot( [coord_exact[0][connect_exact[len(connect)-1][0]] , coord_exact[0][connect_exact[len(connect)-1][1]]] , 
                [value_exact[connect_exact[len(connect)-1][0]] , value_exact[connect_exact[len(connect)-1][1]]] , 
                'b-', linewidth=0.5, label='Exact solution')
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value[connect[len(connect)-1][0]] , value[connect[len(connect)-1][1]]] , 
                'r--', linewidth=0.5 , label='FEM, {} elem'.format(n_elem))
    plt.xlabel('x [m]')
    plt.ylabel('$u_{y0}(x)$ [m]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/u_y0_{}elem.pdf".format(n_elem))
    plt.close()

def PlotTimoshenkoEuler(coord, connect, value_tim, value_eul, L):
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value_tim[connect[i][0]] , value_tim[connect[i][1]]] , 
                'r-', linewidth=0.5)
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value_eul[connect[i][0]] , value_eul[connect[i][1]]] , 
                'b--', linewidth=0.5)
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value_tim[connect[len(connect)-1][0]] , value_tim[connect[len(connect)-1][1]]] , 
                'r-', linewidth=0.5 , label='Timoshenko, {} m'.format(L))
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value_eul[connect[len(connect)-1][0]] , value_eul[connect[len(connect)-1][1]]] , 
                'b--', linewidth=0.5 , label='Euler-Bernoulli, {} m'.format(L))
    plt.xlabel('x [m]')
    plt.ylabel('$u_{y0}(x)$ [m]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/timo_eul_{}.pdf".format(L))
    plt.close()

def PlotTimoshenkoEulerReduced(coord, connect, value_tim, value_eul, L):
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value_tim[connect[i][0]] , value_tim[connect[i][1]]] , 
                'r-', linewidth=0.5)
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value_eul[connect[i][0]] , value_eul[connect[i][1]]] , 
                'b--', linewidth=0.5)
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value_tim[connect[len(connect)-1][0]] , value_tim[connect[len(connect)-1][1]]] , 
                'r-', linewidth=0.5 , label='Timoshenko reduced, {} m'.format(L))
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value_eul[connect[len(connect)-1][0]] , value_eul[connect[len(connect)-1][1]]] , 
                'b--', linewidth=0.5 , label='Euler-Bernoulli, {} m'.format(L))
    plt.xlabel('x [m]')
    plt.ylabel('$u_{y0}(x)$ [m]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/reduced_timo_eul_{}.pdf".format(L))
    plt.close()

def PlotRotation(coord, connect, value, coord_exact, connect_exact, value_exact, n_elem):
    for i in range(len(connect_exact)) : 
        plt.plot( [coord_exact[0][connect_exact[i][0]] , coord_exact[0][connect_exact[i][1]]] , 
                [value_exact[connect_exact[i][0]] , value_exact[connect_exact[i][1]]] , 
                'b-', linewidth=0.5)
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value[connect[i][0]] , value[connect[i][1]]] , 
                'r--', linewidth=0.5)
    plt.plot( [coord_exact[0][connect_exact[len(connect)-1][0]] , coord_exact[0][connect_exact[len(connect)-1][1]]] , 
                [value_exact[connect_exact[len(connect)-1][0]] , value_exact[connect_exact[len(connect)-1][1]]] , 
                'b-', linewidth=0.5, label='Exact solution')
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[len(connect)-1][1]]] , 
                [value[connect[len(connect)-1][0]] , value[connect[len(connect)-1][1]]] , 
                'r--', linewidth=0.5 , label='FEM, {} elem'.format(n_elem))
    plt.xlabel('x [m]')
    plt.ylabel(r'$\theta(x)$ [rad]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/theta_{}elem.pdf".format(n_elem))
    plt.close()

def PlotShearForces(coord, connect, value, coord_exact, connect_exact, value_exact, n_elem):
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value[i][0] , value[i][1]] , 
                'r-', linewidth=0.5 ) 
    for i in range(len(connect_exact)) : 
        plt.plot( [coord_exact[0][connect_exact[i][0]] , coord_exact[0][connect_exact[i][1]]] , 
                [value_exact[connect_exact[i][0]] , value_exact[connect_exact[i][1]]] , 
                'b-', linewidth=0.5 )
    plt.xlabel('x [m]')
    plt.ylabel('V(x) [N]')
    plt.grid(True)
    plt.legend(['FE, {} elem'.format(n_elem),'Exact solution'])
    plt.tight_layout()
    plt.savefig("output/V_{}elem.pdf".format(n_elem))
    plt.close()

def PlotBendingMoments(coord, connect, value, coord_exact, connect_exact, value_exact, n_elem):
    for i in range(len(connect)-1) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [value[i][0] , value[i][1]] , 
                'r-', linewidth=0.5 ) 
    for i in range(len(connect_exact)-1) : 
        plt.plot( [coord_exact[0][connect_exact[i][0]] , coord_exact[0][connect_exact[i][1]]] , 
                [value_exact[connect_exact[i][0]] , value_exact[connect_exact[i][1]]] , 
                'b--', linewidth=0.5 )
    plt.plot( [coord[0][connect[len(connect)-1][0]] , coord[0][connect[connect[len(connect)-1][0]][1]]] , 
                [value[connect[len(connect)-1][0]][0] , value[connect[len(connect)-1][0]][1]] , 
                'r-', linewidth=0.5, label='FE, {} elem'.format(n_elem) )
    plt.plot( [coord_exact[0][connect_exact[len(connect_exact)-1][0]] , coord_exact[0][connect_exact[len(connect_exact)-1][1]]] , 
                [value_exact[connect_exact[len(connect_exact)-1][0]] , value_exact[connect_exact[len(connect_exact)-1][1]]] , 
                'b--', linewidth=0.5 , label='Exact solution')
    plt.xlabel('x [m]')
    plt.ylabel('M(x) [Nm]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/M_{}elem.pdf".format(n_elem))
    plt.close()



# DEFINITIONS
# Display of the shear forces

def PlotShear(coord, connect, Shear) : 
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(coord[0], coord[1], 'ro')
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][0]]+ Shear[i][0], coord[1][connect[i][1]]+ Shear[i][1]] , 
                'g' )
        if not (i == len(connect) or i == 0) : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]]+ Shear[i-1][1], coord[1][connect[i][1]]+ Shear[i][0]] , 
                'g' )
        elif i == 0 : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]], coord[1][connect[i][1]]+ Shear[i][0]] , 
                'g' )
        if i == len(connect) - 1 : 
            plt.plot([coord[0][connect[i][1]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][1]]+ Shear[i][1], coord[1][connect[i][1]]],  
                'g' )
    

def PlotBending(coord, connect, Bending, Shear) : 
    
    for i in range(len(connect)) : 
        plt.plot( [coord[0][connect[i][0]] , coord[0][connect[i][1]]] , 
                [coord[1][connect[i][0]] , coord[1][connect[i][1]]] , 
                '#000000' ) 
    
    plt.plot(coord[0], coord[1], 'ro')
    Starting_point = 0
    
    for i in range(len(connect)) : 
        L = abs(coord[0][connect[i][1]] - coord[0][connect[i][0]])
        x = np.linspace(0, L)
        a =  - (Shear[i][1] - Shear[i][0]) / (2*L)
        b =  - Shear[i][0]
        c = Bending[i][0]
        M = a * x ** 2 + b * x + c
        
        plt.plot(x+Starting_point,M,'r')
        
        Starting_point += L
                     
        if not (i == len(connect) or i == 0) : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]]+ Bending[i-1][1], coord[1][connect[i][1]]+ Bending[i][0]] , 
                'r' )
        elif i == 0 : 
            plt.plot([coord[0][connect[i][0]], coord[0][connect[i][0]]], 
                [coord[1][connect[i][0]], coord[1][connect[i][1]]+ Bending[i][0]] , 
                'r' )
        if i == len(connect) - 1 : 
            plt.plot([coord[0][connect[i][1]], coord[0][connect[i][1]]], 
                [coord[1][connect[i][1]]+ Bending[i][1], coord[1][connect[i][1]]],  
                'r' )
    

"""### 1: INPUT

All parameters that need to be set by the user are listed in this section with # USER
"""

# Parameters of the problem

# Bridge:
E = 30e9 # [N/m^2]
L = 5 # [m]
h = 1 # [m]
b = 0.25 # [m]
A = b*h # [m^2]
I = b * (h**3)/12 # [m^4]
nu = 0.2 #poisson's ratio
k = 5/6 #shear coefficient
n_elem = 4
n_elem_exact = 1000
Timoshenko = True #Timoshenko if true and Bernoulli if false
G = E/(2*(1+nu))  # [N/m^2]
A_c = k*A # [m^2]
reduced = False

# USER
# Coordinates of the nodes: x in the first line and y in the second line, in [m].

# Bridge:
Coord = np.array([np.arange(0,L+1e-10,L/n_elem),np.zeros(n_elem+1)])
Coord_exact = np.array([np.arange(0,L+1e-10,L/n_elem_exact),np.zeros(n_elem_exact+1)])

plt.close(1)
plt.figure(1)
#PlotUndeformed(Coord, np.array([[0,0]])) # Display of nodes

# Total number of degrees of freedom + numbering

No_Ddl = len(Coord[1])*3  # 3 DoF per node

Num_Ddl = np.arange(No_Ddl) # Indexing starts at 0 in Python

print('The structure has ' + str(No_Ddl) + ' degrees of freedom')

# Bridge:
Connect_list = []
for i in range(n_elem):
    Connect_list.append(np.array([i,i+1]))
Connect = np.array(Connect_list)
Connect_list_exact = []
for i in range(n_elem_exact):
    Connect_list_exact.append(np.array([i,i+1]))
Connect_exact = np.array(Connect_list_exact)

#Defining the type of Element : 6 or 5 DoFs (5 - to complete in Assignment)

# Bridge:
Elem_Types = np.ones(n_elem)*6

plt.close(2)
plt.figure(2)          
#PlotUndeformed(Coord, Connect)


# Total number of elements
No_Elem = len(Connect)
No_Elem_exact = len(Connect_exact)
print("The structure is composed of " + str(No_Elem) + " elements.")

# USER
# Fixed degrees of freedom

# Bridge:
Fixed_DoF = np.array([0, 1, 2])

# Free degrees of freedom
Free_DoF = np.delete(Num_Ddl, Fixed_DoF)
        
print(Free_DoF)

# Nodal loads

# Initialization
P = np.zeros(No_Ddl)
P_tim = np.zeros(No_Ddl)  #Timoshenko

# USER
# in [N] and [m]

# Bridge:
P_f = np.zeros(len(Free_DoF))
P_f[len(Free_DoF)-1] = -100e3
P_f[len(Free_DoF)-2] = -20e3
P_f_tim = np.zeros(len(Free_DoF))
P_f_tim[len(Free_DoF)-1] = -100e3
P_f_tim[len(Free_DoF)-2] = -20e3

# Building other vectors:
P[Free_DoF] = P_f
P_tim[Free_DoF] = P_f_tim

# USER

# Bridge:
AE_Elem = np.ones(No_Elem)*A*E
EI_Elem = np.ones(No_Elem)*E*I
GA_c_Elem = np.ones(No_Elem)*G*A_c

# Variable used to magnify the plot of deformed shape, if needed
Scale = 1

"""### 2: COMPUTATIONS

#### Phase 1: Initialization of vectors and matrices
"""

# Initialization of the vector of structural displacements
U = np.zeros(No_Ddl)
U_tim = np.zeros(No_Ddl)  # index 'tim' for Timoshenko
# The displacement vector corresponding to the fixed degrees of freedom is a zero vector
U_d = U[Fixed_DoF]
U_d_tim = U_tim[Fixed_DoF]

# Initialization of the local and global stiffness matrices of the elements, and of the structural stiffness matrix
K_str = np.zeros((No_Ddl, No_Ddl))
k_elem_loc = np.zeros((No_Elem,6,6))  # MODIFIE
k_elem_glob = np.zeros((No_Elem,6,6))
K_str_tim = np.zeros((No_Ddl, No_Ddl))
k_elem_loc_tim = np.zeros((No_Elem,6,6))  # MODIFIE
k_elem_glob_tim = np.zeros((No_Elem,6,6))

# Initialization of the matrices containing : 
# 1. The length of the elements
L_Elem = np.zeros(No_Elem)
# 2. Rotation matrix for each element
r_C = np.zeros((No_Elem,6,6)) #MODIFIE
# 3. Assembly matrix
Assemblage = np.zeros((No_Elem, 6))

"""#### Phase 2: Elements' length, rotation matrices, and assembly matrix
Loop over the elements to calculate their respective lengths, rotation matrices, and assembly matrix
"""

for i in range(No_Elem) : 
    
    # 1. Element's length
    L_x = Coord[0][Connect[i][1]] - Coord[0][Connect[i][0]]  # Length in x
    L_y = Coord[1][Connect[i][1]] - Coord[1][Connect[i][0]]  # Length in y 
    L_Elem[i] = np.sqrt(L_x**2 + L_y**2)
    
     # 2. Rotation matrices [[cos sin 0 0],[0 0 cos sin]]
    sin = L_y / L_Elem[i] # Sine of the rotation angle of truss element i
    cos = L_x / L_Elem[i] # Cosine of the rotation angle of truss element i
    r_C[i] = np.array([[cos, sin, 0, 0, 0, 0],
                       [-sin, cos, 0, 0, 0, 0],
                       [0, 0, 1, 0, 0, 0],
                       [0, 0, 0, cos, sin, 0],
                       [0, 0, 0, -sin, cos, 0],
                       [0, 0, 0, 0, 0, 1]])  #MODIFIE
    
    # Auxiliary matrices for the assembly: positioning of local matrices in the global matrix
    Assemblage[i] = np.array([Connect[i][0]*3, 
                              Connect[i][0]*3+1,
                              Connect[i][0]*3+2,
                              Connect[i][1]*3,
                              Connect[i][1]*3+1,
                              Connect[i][1]*3+2])
    Assemblage = Assemblage.astype(int)

"""#### Phase 3: Computation of local stiffness matrices
Loop through the elements to calculate their respective local stiffness matrices in the local (k_loc) and global (k_glob) reference system, followed by assembly into the structural stiffness matrix 
NOTE: the matrix product is written '@' in Python
"""

for elem in range(No_Elem) : 
    
    # Stiffness matrices in the local reference system, 6 DoF
    if Elem_Types[elem] == 6 and reduced==False: 
        k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, 12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2, 0,  -12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2],
                                     [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem],   0,   -6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem]],
                                     [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, -12*EI_Elem[elem]/L_Elem[elem]**3, -6*EI_Elem[elem]/L_Elem[elem]**2, 0,  12*EI_Elem[elem]/L_Elem[elem]**3,   -6*EI_Elem[elem]/L_Elem[elem]**2],
                                     [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem],  0,    -6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem]]])  #MODIFIE
        #Timoshenko
        k_elem_loc_tim[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0,GA_c_Elem[elem]/L_Elem[elem], GA_c_Elem[elem]/2, 0, -GA_c_Elem[elem]/L_Elem[elem],  GA_c_Elem[elem]/2],
                                     [0,GA_c_Elem[elem]/2, EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GA_c_Elem[elem]/3, 0,   -GA_c_Elem[elem]/2,  -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GA_c_Elem[elem]/6],
                                     [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, -GA_c_Elem[elem]/L_Elem[elem], -GA_c_Elem[elem]/2, 0, GA_c_Elem[elem]/L_Elem[elem], -GA_c_Elem[elem]/2],
                                     [0,GA_c_Elem[elem]/2, -EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GA_c_Elem[elem]/6, 0, -GA_c_Elem[elem]/2,  EI_Elem[elem]/L_Elem[elem]+L_Elem[elem]*GA_c_Elem[elem]/3]])#MODIFIE
    
    
    elif Elem_Types[elem] == 6 and reduced == True:
        k_elem_loc[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, 12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2, 0,  -12*EI_Elem[elem]/L_Elem[elem]**3,  6*EI_Elem[elem]/L_Elem[elem]**2],
                                     [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem],   0,   -6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem]],
                                     [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, -12*EI_Elem[elem]/L_Elem[elem]**3, -6*EI_Elem[elem]/L_Elem[elem]**2, 0,  12*EI_Elem[elem]/L_Elem[elem]**3,   -6*EI_Elem[elem]/L_Elem[elem]**2],
                                     [0, 6*EI_Elem[elem]/L_Elem[elem]**2,   2*EI_Elem[elem]/L_Elem[elem],  0,    -6*EI_Elem[elem]/L_Elem[elem]**2,   4*EI_Elem[elem]/L_Elem[elem]]])  #MODIFIE
        #Timoshenko reduced
        k_elem_loc_tim[elem] = np.array([[AE_Elem[elem]/L_Elem[elem], 0, 0, -AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0,GA_c_Elem[elem]/L_Elem[elem], GA_c_Elem[elem], 0, -GA_c_Elem[elem]/L_Elem[elem],  0],
                                     [0,GA_c_Elem[elem], EI_Elem[elem]/L_Elem[elem]+GA_c_Elem[elem]*L_Elem[elem], 0,   -GA_c_Elem[elem],  -EI_Elem[elem]/L_Elem[elem]],
                                     [-AE_Elem[elem]/L_Elem[elem], 0, 0, AE_Elem[elem]/L_Elem[elem], 0, 0],
                                     [0, -GA_c_Elem[elem]/L_Elem[elem], -GA_c_Elem[elem], 0, GA_c_Elem[elem]/L_Elem[elem], 0],
                                     [0,0, -EI_Elem[elem]/L_Elem[elem], 0, 0,  EI_Elem[elem]/L_Elem[elem]]])#MODIFIE
    

            
    # Stiffness matrices in the global reference system
    k_elem_glob[elem] = np.transpose(r_C[elem]) @ k_elem_loc[elem] @ r_C[elem]
    k_elem_glob_tim[elem] = np.transpose(r_C[elem]) @ k_elem_loc_tim[elem] @ r_C[elem]
    
    # Assembly of the global structural stiffness matrix
    for j in range(len(k_elem_glob[elem])) : 
        for k in range(len(k_elem_glob[elem])) : 
            K_str[Assemblage[elem][j]][Assemblage[elem][k]] += k_elem_glob[elem][j][k]
            K_str_tim[Assemblage[elem][j]][Assemblage[elem][k]] += k_elem_glob_tim[elem][j][k] #Timoshenko

"""#### Phase 4: Partitioning of the stiffness matrix"""

# Sub-matrix for the free DoFs:
K_ff = K_str[Free_DoF[:,None], Free_DoF[None,:]]
K_ff_tim = K_str_tim[Free_DoF[:,None], Free_DoF[None,:]]

#print(K_ff)

# Sub-matrix for the fixed DoFs: 
K_dd = K_str[Fixed_DoF[:,None], Fixed_DoF[None,:]]
K_dd_tim = K_str_tim[Fixed_DoF[:,None], Fixed_DoF[None,:]]

# Sub-matrices K_fd et K_df:
K_fd = K_str[Free_DoF[:,None], Fixed_DoF[None,:]]
K_df = np.transpose(K_fd)
K_fd_tim = K_str_tim[Free_DoF[:,None], Fixed_DoF[None,:]]
K_df_tim = np.transpose(K_fd_tim)

"""#### Phase 5: Displacement's equation

Solving the displacement equation to find the value of the free degrees of freedom
"""

U_f = inv(K_ff) @ (P_f - K_fd @ U_d)
U_f_tim = inv(K_ff_tim) @ (P_f_tim - K_fd_tim @ U_d_tim)

# Completing the global displacement vector:
U[Free_DoF] = U_f
U[Fixed_DoF] = U_d
print(U)
U_tim[Free_DoF] = U_f_tim
U_tim[Fixed_DoF] = U_d_tim

#Computation of exact solution
U_exact = np.zeros(No_Elem+1)
theta_exact = np.zeros(No_Elem+1)
for i in range(No_Elem):
    U_exact[i+1] = -20*(10**3)/(GA_c_Elem[i])*Coord[0][i+1] -100*(10**3)/(EI_Elem[i])*(Coord[0][i+1]**2) + 10*(10**3)/(3*EI_Elem[i])*(Coord[0][i+1]**3)
    theta_exact[i+1] = 10*(10**3)/(EI_Elem[i])*(Coord[0][i+1]**2) - 200*(10**3)/(EI_Elem[i])*Coord[0][i+1]
print(U_exact)
print(theta_exact)

U_exact = np.zeros(No_Elem_exact+1)
theta_exact = np.zeros(No_Elem_exact+1)
V_exact = np.zeros(No_Elem_exact+1)
M_exact = np.zeros(No_Elem_exact+1)
V_exact[0] = -20e3
M_exact[0] = -200e3
for i in range(No_Elem_exact):
    U_exact[i+1] = -20*(10**3)/(GA_c_Elem[0])*Coord_exact[0][i+1] -100*(10**3)/(EI_Elem[0])*(Coord_exact[0][i+1]**2) + 10*(10**3)/(3*EI_Elem[0])*(Coord_exact[0][i+1]**3)
    theta_exact[i+1] = 10*(10**3)/(EI_Elem[0])*(Coord_exact[0][i+1]**2) - 200*(10**3)/(EI_Elem[0])*Coord_exact[0][i+1]
    V_exact[i+1] = -20e3
    M_exact[i+1] = 20*(10**3)*(Coord_exact[0][i+1]-10)

"""#### Phase 6: Reactions' equation

Solving the reactions' equation to find the unknown reactions in the structure
"""

P_d = K_df @ U_f + K_dd @ U_d
P_d_tim = K_df_tim @ U_f_tim + K_dd_tim @ U_d_tim

# Completing the vector of nodal forces:
P[Fixed_DoF] = P_d
P_tim[Fixed_DoF] = P_d_tim

# Computing the structural resisting forces: 
P_r = K_str @ U
P_r_tim = K_str_tim @ U_tim

"""#### Phase 7: Computation of internal forces"""

# Initialization of vectors u_loc et p_loc

u_loc = np.zeros((No_Elem, 6))  #MODIFIE
p_loc = np.zeros((No_Elem, 6))  #MODIFIE
u_loc_tim = np.zeros((No_Elem, 6))  #MODIFIE
p_loc_tim = np.zeros((No_Elem, 6))  #MODIFIE

for i in range(No_Elem) : 
    u_loc[i] = r_C[i] @ U[Assemblage[i]]
    p_loc[i] = k_elem_loc[i] @ u_loc[i]
    u_loc_tim[i] = r_C[i] @ U_tim[Assemblage[i]]
    p_loc_tim[i] = k_elem_loc_tim[i] @ u_loc_tim[i]
    
    #print(p_loc[i])
print('p_loc: {}'.format(p_loc))   
"""### 3: DISPLAY
Display of the internal forces per element
"""

#Print deformed structure
Disp = np.zeros((2,len(Coord[0])))
Disp[0] = U_tim[np.arange(len(Coord[0]))*3]
Disp[1] = U_tim[np.arange(len(Coord[0]))*3+1]
rota = U_tim[np.arange(len(Coord[0]))*3+2]
Disp_Euler = U[np.arange(len(Coord[0]))*3+1]

#PlotDeformed(Coord, Connect, Disp, Scale)
#PlotTransversalDisplacement(Coord, Connect, Disp[1], Coord_exact, Connect_exact, U_exact, n_elem)
#PlotRotation(Coord, Connect, rota, Coord_exact, Connect_exact, theta_exact, n_elem)
"""
if reduced==False:
    PlotTimoshenkoEuler(Coord, Connect, Disp[1], Disp_Euler, L)
if reduced==True:
    PlotTimoshenkoEulerReduced(Coord, Connect, Disp[1], Disp_Euler, L)
"""
#Choice of the element to display 

Elem_ID_to_display = 0

Shear = np.zeros((n_elem,2))
Bending = np.zeros((n_elem,2))
for i in range(n_elem):
    Shear[i][0] = p_loc_tim[i][1]
    Shear[i][1] = -p_loc_tim[i][4]
    Bending[i][0] = p_loc_tim[i][2]
    Bending[i][1] = -p_loc_tim[i][5]
N = p_loc_tim[Elem_ID_to_display][3]

#print(p_loc[Elem_ID_to_display])
#print(Bending)

# Display of axial force
#print("Axial force in element {} = {}kN".format(Elem_ID_to_display, np.around(N/1000,2)))

plt.close(3)
plt.figure(3)
# Display of shear force
PlotShearForces(Coord, Connect, Shear, Coord_exact, Connect_exact, V_exact, n_elem)    
PlotShear(np.array([[0, L_Elem[Elem_ID_to_display]],[0,0]]), np.array([[0,1]]), Shear) 
#print(Shear)

plt.close(4)
plt.figure(4)
# Display of bending moment
PlotBendingMoments(Coord, Connect, Bending, Coord_exact, Connect_exact, M_exact, n_elem) 
PlotBending(np.array([[0, L_Elem[Elem_ID_to_display]],[0,0]]), np.array([[0,1]]), Bending, Shear) 
#print(Bending)





