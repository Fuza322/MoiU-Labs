import numpy as np

def print_result(A_: np.array):
    for j in range(len(A_)):
        print(" ".join((map(str, A_[j]))))

n, i = tuple(map(int, input().split()))

matrix_A = np.zeros((n,n))
matrix_B = np.zeros((n,n))

for j in range(n):
    matrix_A[j] = list(map(float, input().split()))
for j in range(n):
    matrix_B[j] = list(map(float, input().split()))

vect_x = list(map(float, input().split()))
matrix_A[ :, i - 1] = vect_x

vect1_l = matrix_B @ vect_x
if (vect1_l[i - 1] == 0):
    print("NO")
else:
    print("YES")
    vect2_l = vect1_l.copy()
    vect2_l[i - 1] = -1
    vect3_l = vect2_l / (-vect1_l[i - 1])

    matrix_Q = np.identity(len(matrix_A))
    matrix_Q[:, i - 1] = vect3_l[ : ]

    A_ = matrix_Q @ matrix_B

    print_result(A_)