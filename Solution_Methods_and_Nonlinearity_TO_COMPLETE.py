"""LGCIV2041: Course 'Numerical Analysis of Civil Engineering Structures'

Accounting for large nodal displacements - 'Newton-Raphson method' and 'Displacement-control method' for column with nonlinear spring at the base
March 2023, João Pacheco de Almeida, Hadrien Rattez, Alexandre Sac-Morane

The 2 parts of the code that should be completed are indicated below by:
================================= PART 1 TO COMPLETE: START HERE =================================================
================================= PART 1 TO COMPLETE: FINISH HERE ===================================
================================= PART 2 TO COMPLETE: START HERE =================================================
================================= PART 2 TO COMPLETE: FINISH HERE ==================================="""

import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin,pi, atan2

def Newton_Raphson_Ramberg_Osgood(theta_applied, tolerance_theta, theta_y, gamma, M_y):
    
    # Function that performs a Newton-Rhapson method to obtain the resisting moment M and the stiffness k_PH of a plastic hinge following the Ramberg-Osgood law

    # Convergence flag initialization (false => not converged; true => converged)
    conv = False
    # Initialization of the iteration counter:
    iteration = 0
    # Initialization of moment for iterative procedure:
    M = 0
    # The Ramberg-Osgood expression is intended to be used for positive values of the moment-theta relation. 
    # Hence, when a negative value of theta_applied, the solver finds with the NR method the values of moment + stiffness for
    # the symmetric value of theta_applied (i.e., positive), and after it changes the sign of such moment (the stiffness is the same)
    theta_applied_negative = False
    if theta_applied < 0:
        theta_applied_negative = True
        theta_applied = -theta_applied

    while conv == False: # Newton-Raphson iterations:
        # Compute theta:
        theta = theta_y*(M/M_y*(1+(M/M_y)**(gamma-1)))  
        # Compute stiffness:
        dtheta_dM = (theta_y*((M/M_y)**(gamma - 1) + 1))/M_y + (M*theta_y*(M/M_y)**(gamma - 2)*(gamma - 1))/M_y**2
        k_PH = 1 / dtheta_dM
        # Compute residual:
        Res = theta_applied - theta # Residual
        #  Compute residual norm for convergence:
        Residual = np.linalg.norm(Res) # Euclidean norm of residual vector
        # Check for convergence:
        if Residual <= tolerance_theta: # Check for convergence
            conv = True # Iterative process has converged
        else:
            M = M + Res / dtheta_dM
        # Update iteration:
        iteration = iteration + 1
        
    # Adjustment for the sign of M as mentioned above...:
    if theta_applied_negative == True:
        M = -M

    return M, k_PH

# ====================================
# ====== INPUT PARAMETERS: START ======
# ====================================
# Account for nonlinear geometric effects (large displacements):
Consider_non_linear_geometric_effects = False
# Use classical Newton-Raphson (1) or Displacement-Control method (2):
Classical_NR_or_Disp_Control = 2
# Increments of lateral displacement in [m] (for Newton-Raphson method)
F_lat = np.linspace(0, 300.0, num = 46) # Use for NR with linear geometry
# Increments of lateral displacement in [m] (for Displacement-Control method)
Delta_lat = np.linspace(0.0, -0.8, num=24)
# Max no. of iterations for solver
Max_no_iterations = 15    

# Input data:
EI = 2*10**8 # kNm2
EA = 1.5*10**7  # kN
L = 5 # m
W = 1000 # kN
theta_y = 0.035 # rad
gamma = 5
M_y = 3000  # kNm
tol_force = 1 # kN
tol_M_PH = 1*10**-5 # rad tolerance for the Ramberg-Osgood Newton-Raphson procedure
# ==================================
# ====== INPUT PARAMETERS: END ======
# ===================================

# Setting dimensions of variables:
if Classical_NR_or_Disp_Control == 1: # Using classical Newton-Raphson method:
    No_increments = len(F_lat)
else:  # Using displacement-control method:
    No_increments = len(Delta_lat)    

U_conv = np.zeros((4, No_increments),dtype = float)

# Assigning the imposed displacement to the second DoF:
if Classical_NR_or_Disp_Control == 2:
    U_conv[1,] = Delta_lat

P_r_conv = np.zeros((4, No_increments), dtype = float)
K_str_conv = np.zeros((4, 4, No_increments), dtype = float)
Counter_Iterations = np.zeros((2, No_increments), dtype = float) # variable that stores the number of iterations per increment

# Next commands are really not necessary, just filling in a first line with the imposed force / displacements
# Further down the code we will add a 2nd line with the number of iterations for each of those imposed force / displacements
if Classical_NR_or_Disp_Control == 1: # Using classical Newton-Raphson method:
    Counter_Iterations[0, ] = F_lat
else:  # Using displacement-control method:
    Counter_Iterations[0, ] = Delta_lat

# Initializing residual:
Res = np.zeros((4),dtype = float)

# Initial angle of the column:
Theta_column = pi/2
c=cos(Theta_column)
s=sin(Theta_column)
# Compute compatibility matrix between the global-local reference system:
Compatibility_matrix_global_local = np.array([[c, s, 0, 0, 0, 0], 
                                              [-s,c, 0, 0, 0, 0],
                                              [0, 0, 1, 0, 0, 0],
                                              [0, 0, 0, c, s, 0],
                                              [0, 0, 0,-s, c, 0],
                                              [0, 0, 0, 0, 0, 1]])
# Compute equilibrium matrix between the local-global reference system:
Equilibrium_matrix_local_global = np.transpose(Compatibility_matrix_global_local)

## Cycle through the load increments:
for i in range(No_increments):
    
    # Initialization of displacement vector for iterative procedure:
    if i ==0:
        U = U_conv[:, i].copy()
    else:
        U = U_conv[:, i-1].copy()
        if Classical_NR_or_Disp_Control == 2:
            U[1]=U_conv[1, i]
    # Convergence flag initialization (false => not converged true => converged):
    conv = False   
    # Initialization of the iteration counter:
    iteration = 0     
    
    ## global NR procedure
    while (conv == False and (iteration <= Max_no_iterations)): # Newton-Raphson iterations
        
        # State Determination - Computation of structural resisting forces (in the global reference system):
        
        if Consider_non_linear_geometric_effects == False:  # If linear geometry is considered:
            # ==================================================================================
            #  ===================== COMPUTE P_r with LINEAR GEOMETRY - START =====================
            # ==================================================================================
            # Compute nodal element displacements in the global reference system from the nodal structural displacements (in the global reference system):
            u_global = np.hstack((np.array([0, 0]), U))
            # Compute nodal element displacements in the local reference system:
            u_loc = np.dot(Compatibility_matrix_global_local, u_global)
            # Compute basic element displacements (in the basic reference system) considering linear geometry:
            Compatibility_matrix_local_basic_linear = np.array([[ -1,  0,  0, 1,   0,  0],
                                                                [  0, 1/L, 1, 0, -1/L, 0],
                                                                [  0, 1/L, 0, 0, -1/L, 1]])
            u_bsc = np.dot(Compatibility_matrix_local_basic_linear, u_loc)
            # Compute basic element forces from basic displacements, using the basic stiffness matrix k_bsc:
            k_bsc = np.array([[ EA/L ,  0   ,  0   ],
                              [   0   ,4*EI/L,2*EI/L],
                              [   0   ,2*EI/L,4*EI/L]])
            p_bsc = np.dot(k_bsc,  u_bsc)
            # Compute nodal element forces in the local reference system from basic element forces (in the basic reference system):
            Equilibrium_matrix_basic_local_linear = np.transpose(Compatibility_matrix_local_basic_linear)
            p_loc = np.dot(Equilibrium_matrix_basic_local_linear, p_bsc)
            # Compute nodal element forces in the global reference system:
            p_global = np.dot(Equilibrium_matrix_local_global, p_loc)
            # Compute the structural resisting forces P_r (in the global reference system), which are composed of an elastic contribution from the element and another from the nonlinear spring:
            P_r_elastic = p_global[2:] # Compute P_r_elastic (contribution from the elastic element)
            M_PH, k_PH = Newton_Raphson_Ramberg_Osgood(U[0], tol_M_PH, theta_y, gamma, M_y) # Calling Newton-Rhapson method to find the moment and stiffness corresponding to a rotation U(1) from the Ramberg-Osgood formula (contribution from the nonlinear spring)
            P_r = P_r_elastic + np.hstack((M_PH, np.array([0, 0, 0]))) # Compute P_r (sum of the elastic contribution from the element and the nonlinear spring)
            # ==================================================================================            
            #  ===================== COMPUTE P_r with LINEAR GEOMETRY - END =====================
            # ==================================================================================            
        else :        # If nonlinear geometric effects (large displacements) are considered:
            # ==================================================================================================
            # ================================= PART 1 TO COMPLETE: START HERE ===================================
            #  ========================== COMPUTE P_r with NONLINEAR GEOMETRIC EFFECTS - START ====================
            # ===================================================================================================       
            # Introduce below the commands that allow to obtain the basic element displacements (in the basic reference system) from the nodal structural displacements (in the global reference system), considering nonlinear geometry:

            
            
            
            
            # Compute basic forces from basic displacements, using the basic stiffness matrix k_bsc:





            # Introduce below the commands that allow to compute the nodal element forces in the global reference system from the basic element forces:





            # Compute the structural resisting forces P_r (in the global reference system), which are composed of an elastic contribution from the element and another from the nonlinear spring:
            P_r_elastic = p_global[2:] # Compute P_r_elastic (contribution from the elastic element)
            M_PH, k_PH = Newton_Raphson_Ramberg_Osgood(U[0], tol_M_PH, theta_y, gamma, M_y) # Calling Newton-Rhapson method to find the moment and stiffness corresponding to a rotation U(1) from the Ramberg-Osgood formula (contribution from the nonlinear spring)
            P_r = P_r_elastic + np.hstack((M_PH, np.array([0, 0, 0]))) # Compute P_r (sum of the elastic contribution from the element and the nonlinear spring)         
            # ==================================================================================================
            # ================================= PART 1 TO COMPLETE: FINISH HERE ===================================
            #  ========================== COMPUTE P_r with NONLINEAR GEOMETRIC EFFECTS - END ======================            
            # ===================================================================================================      
        
        
        # State Determination - Computation of the structural stiffness matrix in the global reference system:
        
        if Consider_non_linear_geometric_effects == False:         # If linear geometry is considered:
            # ==================================================================================
            #  ===================== COMPUTE k_loc with LINEAR GEOMETRY - START =====================
            # ==================================================================================
            # Compute element stiffness matrix in the local reference system:
            k_loc = np.dot(Equilibrium_matrix_basic_local_linear, np.dot( k_bsc, Compatibility_matrix_local_basic_linear))
            
            # Computing the global element stiffness from the local element stiffness matrix:
            k_global = np.dot(Equilibrium_matrix_local_global, np.dot(k_loc, Compatibility_matrix_global_local))
            # The plastic hinge stiffness k_PH was computed above...
            # Structural stiffness:
            K_str = k_global[2:,2:]
            K_str[0,0] = K_str[0,0] + k_PH # Adding the stiffness from the nonlinear spring...
            # ==================================================================================
            #  ===================== COMPUTE k_loc with LINEAR GEOMETRY - END =====================
            # ==================================================================================
        else:             # If nonlinear geometric effects (large displacements) are considered:
            # ==================================================================================================
            # ================================= PART 2 TO COMPLETE: START HERE ===================================
            #  ========================== COMPUTE k_loc with NONLINEAR GEOMETRIC EFFECTS - START ====================
            # ===================================================================================================                
            # Compute local stiffness matrix:
            
            

            # Computing the global element stiffness from the local element stiffness matrix:

            
            
            # The plastic hinge stiffness k_PH was computed above...
            # Structural stiffness:
            K_str = k_global[2:,2:]
            K_str[0,0] = K_str[0,0] + k_PH # Adding the stiffness from the nonlinear spring...
            # ==================================================================================================
            # ================================= PART 2 TO COMPLETE: FINISH HERE ===================================
            #  ========================== COMPUTE k_loc with NONLINEAR GEOMETRIC EFFECTS - END ======================            
            # ===================================================================================================                

        # Evaluate convergence + Solve linearized system (displacement control) + Update displacements:

        if Classical_NR_or_Disp_Control == 1: # Use classical Newton-Raphson method
            # Compute residual:
            P = [0, -F_lat[i], -W, 0]
            Res = P - P_r # Residual
            #  Compute residual norm for convergence:
            Residual = np.linalg.norm(Res) # Euclidean norm of residual vector
            if Residual <= tol_force: # Check for convergence
                conv = True # Iterative process has converged
            else:
                U = U + np.linalg.solve(K_str, Res)
            # Check for convergence
        else:  # Use displacement-control method
            # Compute residual:
            P = np.array([0, -W, 0])

            aux = np.array([0,2,3])
            Res[np.ix_(aux)] = P- P_r[np.ix_(aux)]
            # Residual for force-controlled component
            Res[1] = 0 # Residual for the displacement-controlled component
            #  Compute residual norm for convergence:
            Residual = np.linalg.norm(Res) # Euclidean norm of residual vector
            if Residual <= tol_force: # Check for convergence
                conv = True # Iterative process has converged
            else:
                K_str[1, :] = 0 # zeroing columns corresponding to displacement-controlled dof
                K_str[:, 1] = 0 # zeroing rows corresponding to displacement-controlled dof
                K_str[1, 1] = 1 # Assigning ones at diagonal entries corresponding to displacement-controlled dof
                U = U + np.linalg.solve(K_str, Res)
                # Check for convergence
       
        # Update iteration:
        iteration = iteration + 1
        
    Counter_Iterations[1, i] = iteration
    U_conv[:, i] = U
    P_r_conv[:, i] = P_r
    K_str_conv[:, :, i] = K_str


## ---- Plotting ----
fig = plt.figure()  
plt.plot(-U_conv[1,:],-P_r_conv[1,:], 'ro-')
np.savetxt('Displacement.out', (-np.transpose(U_conv[1,:])))   # x,y,z equal sized 1D arrays
np.savetxt('Force.out', (np.transpose(-P_r_conv[1,:])))
plt.title('Force-displacement response')
plt.xlabel('Lateral Displacement [m]')
plt.ylabel('Lateral Force [kN]')

