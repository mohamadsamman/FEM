import numpy as np

# -----------------------------------------------------
# Functions for Base Functions and Gradients (P1 elements)
# -----------------------------------------------------

def Lambda_of_P1(x, x1, x2):
    Lambda1 = (x - x2) / (x1 - x2)
    Lambda2 = (x - x1) / (x2 - x1)
    return np.array([Lambda1, Lambda2])

def MGrad_Lambda_of_P1(x1, x2):
    G1 = 1.0 / (x1 - x2)
    G2 = 1.0 / (x2 - x1)
    return np.array([G1, G2])

# Function for validation of numerical integration
def MyF(x):
    return x

def MyA11(L1, L2, dl1, dl2):
    return dl1 * dl1 + L1 * L1

def MyA12(L1, L2, dl1, dl2):
    return dl1 * dl2 + L1 * L2

def MyA21(L1, L2, dl1, dl2):
    return dl2 * dl1 + L2 * L1

def MyA22(L1, L2, dl1, dl2):
    return dl2 * dl2 + L2 * L2

def MyF1(L1, L2, x):
    return L1 * MyF(x)

def MyF2(L1, L2, x):
    return L2 * MyF(x)

# Integration Functions (Numerical Integration of Functions)
def IntegElmtA(x1, x2, Ng, Omega, Lambda1, Lambda2, dl1, dl2, MyFm):
    out = 0
    for g in range(Ng):
        out += Omega[g] * MyFm(Lambda1[g], Lambda2[g], dl1, dl2)
    out *= abs(x1 - x2)
    return out

def IntegElmtF(x1, x2, Ng, Omega, Lambda1, Lambda2, MyFm):
    out = 0
    for g in range(Ng):
        xg = x1 * Lambda1[g] + x2 * Lambda2[g]
        out += Omega[g] * MyFm(Lambda1[g], Lambda2[g], xg)
    out *= abs(x1 - x2)
    return out

# ----------------------------------------------------------------------
# Gauss Quadrature Points and Weights for Numerical Integration
# ----------------------------------------------------------------------
def gauss_quadrature(Ng):
    if Ng == 1:
        Omega = np.array([1.0])
        Lambda = np.array([[0.5]])
    elif Ng == 2:
        Omega = np.array([0.5, 0.5])
        Lambda = np.array([[0.5 * (1 - 1/np.sqrt(3))],
                           [0.5 * (1 + 1/np.sqrt(3))]])
    elif Ng == 3:
        Omega = np.array([5.0/18.0, 8.0/18.0, 5.0/18.0])
        Lambda = np.array([[0.5 * (1 - np.sqrt(3/5))],
                           [0.5],
                           [0.5 * (1 + np.sqrt(3/5))]])
    Lambda2 = 1 - Lambda
    return Omega, Lambda, Lambda2

# Main Execution for Assembly and Solution

def assemble_element(x1e, x2e, Ng, Omega, Lambda1, Lambda2, dl1dx, dl2dx):
    AeP1 = np.zeros((2, 2))
    AeP1[0, 0] = IntegElmtA(x1e, x2e, Ng, Omega, Lambda1, Lambda2, dl1dx, dl2dx, MyA11)
    AeP1[0, 1] = IntegElmtA(x1e, x2e, Ng, Omega, Lambda1, Lambda2, dl1dx, dl2dx, MyA12)
    AeP1[1, 0] = IntegElmtA(x1e, x2e, Ng, Omega, Lambda1, Lambda2, dl1dx, dl2dx, MyA21)
    AeP1[1, 1] = IntegElmtA(x1e, x2e, Ng, Omega, Lambda1, Lambda2, dl1dx, dl2dx, MyA22)
    
    BeP1 = np.zeros(2)
    BeP1[0] = IntegElmtF(x1e, x2e, Ng, Omega, Lambda1, Lambda2, MyF1)
    BeP1[1] = IntegElmtF(x1e, x2e, Ng, Omega, Lambda1, Lambda2, MyF2)
    
    return AeP1, BeP1

# Define parameters
x1e = 2  # Example value
x2e = 3  # Example value
Ng = 3  # Number of Gauss points

# Compute Gauss points and weights
Omega, Lambda, Lambda2 = gauss_quadrature(Ng)

# Compute gradients for basis functions
Grad_Lambda = MGrad_Lambda_of_P1(x1e, x2e)
dl1dx = Grad_Lambda[0]
dl2dx = Grad_Lambda[1]

# Assemble element stiffness matrix and force vector
AeP1, BeP1 = assemble_element(x1e, x2e, Ng, Omega, Lambda[0], Lambda2[0], dl1dx, dl2dx)

# Analytical solution for validation
Aex = np.array([[4/3, -5/6], [-5/6, 4/3]])
Bex = np.array([7/6, 8/6])

# Display norms of the errors
print(f"Error in A: {np.linalg.norm(AeP1 - Aex)}")
print(f"Error in B: {np.linalg.norm(BeP1 - Bex)}")

# Boundary conditions and solution (optional)
Ibord = 2  # Example boundary condition selection
if Ibord == 1:
    Exact = np.array([0, 0, 3/16])
elif Ibord == 2:
    Exact = np.array([2/3, 0, 25/48])

# Solve the system (example)
Mat_Ae = np.array([[4/3, -5/6], [-5/6, 4/3]])  # Example matrix
Vec_Be = np.array([7/6, 8/6])  # Example vector
Resu = np.linalg.solve(Mat_Ae, Vec_Be)  # Solve Ax = B
print(f"Error in solution: {np.linalg.norm(Resu - Exact)}")

