import math
import numpy as np
 
def simplex_method_Main_Phase(vect_c: np.array, A: np.array, A_rev_basic: np.array, vect_x: np.array, j: np.array):
    A_basic = np.zeros((M, len(j)))
    
    for k in range(len(J)):
        A_basic[: , k] = A[: , j[k] - 1]
    
    if A_rev_basic is None:
        A_rev_basic = np.linalg.inv(A_basic)
    while True:
        
        cb = np.array([vect_c[i - 1] for i in j])
        
        vect_u = np.array(cb @ A_rev_basic)
        
        Δ = np.array(vect_u @ A) - vect_c
        
        J0 = -1    
        min_Δ = math.inf
        for i in range(len(vect_c)):
            if i + 1 not in J:
                if Δ[i] < 0 and Δ[i] < min_Δ:
                    J0 = i
                    min_Δ = Δ[i]
        if J0 == -1:
            print("Bounded")
            for i in range(len(vect_x)):
                print(f"{vect_x[i]} ")
            return
        vect_z = A_rev_basic @ A[: , J0]
        
        tetes = [float(vect_x[J[i] - 1] / vect_z[i]) if vect_z[i] > 0 else math.inf for i in range(len(vect_z))]
        min_teta = min(tetes)
        if min_teta == math.inf:
            print("Unbounded")
            return
        
        min_J = tetes.index(min_teta)
        J[min_J] = J0 + 1
        
        updated_x = np.zeros(len(vect_x), dtype = float)
        updated_x[J - 1] = vect_x[j - 1] - min_teta * vect_z
        updated_x[J0] = min_teta
        vect_x = updated_x
        
        A_rev_basic = sherman_morrison(A_rev_basic, A_basic, A[: , J0], min_J + 1)
        if A_rev_basic is None:
            print("Unbounded")
            return

def sherman_morrison(matrix_B: np.array, A: np.array, vect_x: np.array, i: int) -> np.array:
    A_ = A.copy()
    A_[: , i - 1] = vect_x
 
    vect1_l = matrix_B @ vect_x
    if (vect1_l[i - 1] == 0):
        return None
    else:
        vect2_l = vect1_l.copy()
        vect2_l[i - 1] = -1
        vect3_l = vect2_l / (-vect1_l[i - 1])
        matrix_Q = np.identity(len(A_))
        matrix_Q[: , i - 1] = vect3_l[:]
        A_rev = matrix_Q @ matrix_B
        return A_rev
#----------------------------------------------------------------
N, M = tuple(map(int, input().split()))
A = np.zeros((N,M))
 
for i in range(M):
    A[i] = np.array([*map(float, input().split())], dtype = float)
B = np.array([*map(float, input().split())], dtype = float)
vect_c = np.array([*map(float, input().split())], dtype = float)
vect_x = np.array([*map(float, input().split())], dtype = float)
J = np.array([*map(int, input().split())])
simplex_method_Main_Phase(vect_c, A, None, vect_x, J)