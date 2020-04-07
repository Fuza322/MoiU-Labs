using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace FirstPhase
{
    public static class Program
    {   // Enter the initial conditions
        static double[] InputDouble() =>
           Console.ReadLine().Split(' ').Select(x => double.Parse(x, CultureInfo.InvariantCulture)).ToArray();
        
        static int[] InputInt() =>
            Console.ReadLine().Split(' ').Select(x => int.Parse(x)).ToArray();
     
        static double[,] InputMatrix(int n, int m)
        {
            var matrix = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                var input = InputDouble();

                for (int j = 0; j < m; j++)
                {
                    matrix[i, j] = input[j];
                }
            }
            return matrix;
        }

        // Geting a column of matrix
        static double[] Get_Column(double[,] matrix, int i)
        {
            var n = matrix.LengthRows();
            var result = new double[n];
            for (int j = 0; j < n; j++)
            {
                result[j] = matrix[j, i];
            }
            return result;
        }

        // Geting a row of matrix
        static double[] Get_Row(double[,] matrix, int i)
        {
            var n = matrix.LengthColumns();
            var result = new double[n];
            for (int j = 0; j < n; j++)
            {
                result[j] = matrix[i, j];
            }
            return result;
        }

        // Rows count calculation
        public static int LengthRows(this double[,] matrix)
        {
            return matrix.GetUpperBound(0) + 1;
        }

        // Columns count calculation
        public static int LengthColumns(this double[,] matrix)
        {
            return matrix.GetUpperBound(1) + 1;
        }

        // Multiply the vector by the matrix
        static double[] MultMatrixOnVector(double[,] matrix, double[] vector)
        {
            var result = new double[matrix.LengthRows()];
            for (int i = 0; i < matrix.LengthRows(); i++)
            {
                result[i] = 0;
                for (int j = 0; j < matrix.LengthColumns(); j++)
                {
                    result[i] += matrix[i, j] * vector[j];
                }
            }
            return result;
        }

        // Multiply the matrix by the vector
        static double[] MulMatrixOnVector(double[] vector, double[,] matrix)
        {
            var result = new double[matrix.LengthColumns()];
            for (int i = 0; i < matrix.LengthColumns(); i++)
            {
                result[i] = 0;
                for (int j = 0; j < matrix.LengthRows(); j++)
                {
                    result[i] += matrix[j, i] * vector[j];
                }
            }
            return result;
        }

        // Vector multiplication
        static double MulVectors(double[] vectorA, double[] vectorB)
        {
            double result = 0;
            for (int i = 0; i < vectorA.Length; i++)
            {
                result += vectorA[i] * vectorB[i];
            }
            return result;
        }

        // Sherman-Morrison's formul 
        public static double[,] Sherman_morrison(double[,] matrix, double[] column, int colIndex)
        {
            var n = matrix.LengthRows();
            var vect1_l = MultMatrixOnVector(matrix, column);
            var vect2_l = new double[n];
            vect1_l.CopyTo(vect2_l, 0);
            vect2_l[colIndex] = -1;
            var matrix_Q = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                matrix_Q[i, i] = 1;
                vect2_l[i] = -1 / vect1_l[colIndex] * vect2_l[i];
                matrix_Q[i, colIndex] = vect2_l[i];
            }
            //------------------------------------------------
            var matrix_C = new double[matrix_Q.LengthRows(), matrix.LengthColumns()];
            for (var i = 0; i < colIndex; i++)
            {
                for (var j = 0; j < matrix.LengthColumns(); j++)
                {
                    matrix_C[i, j] = matrix_Q[i, colIndex] * matrix[colIndex, j] + matrix_Q[i, i] * matrix[i, j];
                }
            }
            for (var j = 0; j < matrix.LengthColumns(); j++)
            {
                matrix_C[colIndex, j] = matrix_Q[colIndex, colIndex] * matrix[colIndex, j];
            }
            for (var i = colIndex + 1; i < matrix_Q.LengthRows(); i++)
            {
                for (var j = 0; j < matrix.LengthColumns(); j++)
                {
                    matrix_C[i, j] = matrix_Q[i, colIndex] * matrix[colIndex, j] + matrix_Q[i, i] * matrix[i, j];
                }
            }
            return matrix_C;
        }

        // Transpose matrix
        public static void InverseMatrix(double[,] matrix)
        {
            var n = matrix.LengthRows();
            for (int i = 1; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    var tmp = matrix[i, j];
                    matrix[i, j] = matrix[j, i];
                    matrix[j, i] = tmp;
                }
            }
        }

        // The formation of the basis matrix
        static double[,] CalculateBasisMatrix(double[,] matrix, int[] basis)
        {
            int n = basis.Length;
            var basisMatrix = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    basisMatrix[i, j] = matrix[i, basis[j]];
                }
            }
            return basisMatrix;
        }

        // Getting the basis vector
        static double[] CaltulateBasisVector(double[] vector, int[] basis) =>
            basis.Select(e => vector[e]).ToArray();

        //// Getting the nobasis vector
        static int[] CalculateNonBasis(int[] basis, int n)
        {
            var lBasis = basis.ToHashSet();
            var notBasisVector = new List<int>();
            for (int i = 0; i < n; i++)
            {
                if (!lBasis.Contains(i))
                {
                    notBasisVector.Add(i);
                }
            }
            return notBasisVector.ToArray();
        }

        // Geting the Delta
        static double[] GetDelta(double[] u, double[,] A, double[] c, int[] notBasis)
        {
            var result = new double[A.LengthColumns()];
            foreach (var i in notBasis)
            {
                result[i] = 0;
                for (int j = 0; j < A.LengthRows(); j++)
                {
                    result[i] += A[j, i] * u[j];
                }
                result[i] -= c[i];
            }
            return result;
        }
       
        // Deltas rating
        static int? CheckDelta(double[] Δ, int[] notBasis)
        {
            foreach (var i in notBasis)
            {
                if (Δ[i] < 0)
                { return i; }
            }
            return null;
        }

        // Getting the Z-vector
        static double[] GetVect_Z(double[,] B, double[,] A, int badIndex)
        {
            var result = new double[B.LengthRows()];
            for (int i = 0; i < B.LengthRows(); i++)
            {
                result[i] = 0;
                for (int j = 0; j < B.LengthColumns(); j++)
                {
                    result[i] += B[i, j] * A[j, badIndex];
                }
            }
            return result;
        }

        // Finding Teta 
        static double[] GetTeta(double[] x, int[] basis, double[] z)
        {
            var teta = new double[z.Length];
            for (int i = 0; i < z.Length; i++)
            {
                teta[i] = z[i] > 0 ? x[basis[i]] / z[i] : double.PositiveInfinity;
            }

            return teta;
        }

        // Get minimal Teta with his Index
        static (double, int) GetTeta_min(double[] teta)
        {
            var min = teta[0];
            int minIndex = 0;
            for (int i = 1; i < teta.Length; i++)
            {
                if (teta[i] < min)
                {
                    min = teta[i];
                    minIndex = i;
                }
            }
            return (min, minIndex);
        }

        // Finding a new basis plan
        static (double[], int[]) FindNewPlan(double[] oldPlan, int[] oldBasis, double min_teta, int minIndex, int J_0, double[] z)
        {
            var n = oldPlan.Length;
            var m = oldBasis.Length;
            var newPlan = new double[n];

            newPlan[J_0] = min_teta;
            for (int i = 0; i < m; i++)
            {
                newPlan[oldBasis[i]] = oldPlan[oldBasis[i]] - min_teta * z[i];
            }
            var newBasis = new int[m];
            oldBasis.CopyTo(newBasis, 0);
            newBasis[minIndex] = J_0;
            return (newPlan, newBasis);
        }

        // Algorithm of the main phase of the simplex-method
        static (double[], int[], double[,]) MainPhase(double[,] matrix_A, double[] vect_c, double[] vect_x, int[] A_basis, int n) //
        {
            double min_teta;
            int min_J = 0;
            double[,] A_rev_basis = null;
            int? J_0 = 0;
            while (true)
            {
                var noN_Basis = CalculateNonBasis(A_basis, n);
                double[] vect_u = null;
                if (A_rev_basis == null)
                {
                    A_rev_basis = CalculateBasisMatrix(matrix_A, A_basis);
                    vect_u = CaltulateBasisVector(vect_c, A_basis);
                }
                else
                {
                    A_rev_basis = Sherman_morrison(A_rev_basis, Get_Column(matrix_A, J_0.Value), min_J);
                    vect_u = MulMatrixOnVector(CaltulateBasisVector(vect_c, A_basis), A_rev_basis);
                }
                var Δ = GetDelta(vect_u, matrix_A, vect_c, noN_Basis);
                J_0 = CheckDelta(Δ, noN_Basis);
                if (J_0 == null)
                {
                    return (vect_x, A_basis, A_rev_basis);
                }
                var vect_z = GetVect_Z(A_rev_basis, matrix_A, J_0.Value);
                if (vect_z.All(p => p < 1E-6))
                {
                    return (null, null, null);
                }
                var teta = GetTeta(vect_x, A_basis, vect_z);
                (min_teta, min_J) = GetTeta_min(teta);
                (vect_x, A_basis) = FindNewPlan(vect_x, A_basis, min_teta, min_J, J_0.Value, vect_z);
            }
        }

        // Geting the new matrix after removing the K-row and K-column from the initial solutions
        public static double[,] GetNewMatrix(double[,] matrix, int rowToDelete, int columnForDelete)
        {
            var newMatrix = new double[matrix.LengthRows() - 1, matrix.LengthColumns() - 1];
            for (int i = 0; i < newMatrix.LengthRows(); i++)
            {
                var rowIndex = i >= rowToDelete ? i + 1 : i;
                if (rowIndex == matrix.LengthRows())
                {
                    break;
                }
                for (int j = 0; i < newMatrix.LengthColumns(); j++)
                {
                    var columnIndex = j >= columnForDelete ? j + 1 : j;
                    if (columnIndex == matrix.LengthColumns())
                    {
                        break;
                    }
                    newMatrix[i, j] = matrix[rowIndex, columnIndex];
                }
            }
            return newMatrix;
        }

        static void Main()
        {
            var input = InputInt();

            var m = input[0];
            var n = input[1];
            var A = InputMatrix(m, n);
            var b = InputDouble();
            var c = InputDouble();

            var matrixExpanded = new double[m, n + m];

            c = new double[n + m];

            for (int i = 0; i < m; i++)
            {
                if (b[i] < 0)
                {
                    b[i] = -b[i];
                    for (int j = 0; j < n; j++)
                    {
                        A[i, j] = -A[i, j];
                    }
                }
                for (int j = 0; j < n; j++)
                {
                    matrixExpanded[i, j] = A[i, j];
                }
            }
            var x = new double[n + m];
            var basis = new int[m];
            
            for (int i = n; i < n + m; i++)
            {
                c[i] = -1;
                matrixExpanded[i - n, i] = 1;
                x[i] = b[i - n];
                basis[i - n] = i;
            }
            double[,] AbInv = null;
            
            (x, basis, AbInv) = MainPhase(matrixExpanded, c, x, basis, n + m);

            if (x == null)
            {
                Console.WriteLine("Unbounded");
                return;
            }
            // 
            for (int i = n; i < n + m; i++)
            {
                if (Math.Abs(x[i]) > 1E-6)
                {
                    Console.WriteLine("No solution");
                    return;
                }
            }

            var replaced = false;
            var isFirstCheck = true;
            int k = 0;
            double[] newColumn = null;
            while (true)
            {
                //
                if (basis.All(e => e < n))
                {
                    Console.WriteLine("Bounded");
                    //x = x[..n];
                    x = x.Take(n).ToArray();
                    for (var i = 0; i < x.Length; i++)
                    {
                        Console.Write(x[i].ToString(CultureInfo.GetCultureInfo("en-US")) + " ");
                    }
                    return;
                }

                //
                var Ab = CalculateBasisMatrix(matrixExpanded, basis);
                if (!isFirstCheck)
                {
                    if (!replaced)
                    {
                        // Finding a new inversion after we throw out the row and column
                        var e = new double[AbInv.LengthRows()];
                        e[k] = 1;
                        var matrix = Sherman_morrison(AbInv, e, k); 
                        InverseMatrix(matrix);
                        matrix = Sherman_morrison(AbInv, e, k); 
                        var newMatrix = GetNewMatrix(matrix, k, k);
                        AbInv = newMatrix;
                    }
                    else
                    {
                        AbInv = Sherman_morrison(AbInv, newColumn, k);
                    }
                }
                isFirstCheck = false;

                var notBasis = CalculateNonBasis(basis, n);
                k = 0;
                
                for (int i = 0; i < m; i++)
                {
                    if (basis[i] > n)
                    {
                        k = i; 
                        break;
                    }
                }

                if (!replaced)
                {
                    matrixExpanded = GetNewMatrix(matrixExpanded, basis[k] - n, basis[k]);
                    var newBasis = basis.ToList();
                    for (int i = k + 1; i < m; i++)
                    {
                        newBasis[i]--;
                    }
                    m--;
                    newBasis.RemoveAt(k);
                    basis = newBasis.ToArray();
                }

                replaced = false;
                foreach (var i in notBasis)
                {
                    if (i < Ab.LengthColumns())
                    {
                        var col = Get_Column(Ab, i);
                        var alpha = MulVectors(Get_Row(AbInv, k), col);
                        if (Math.Abs(alpha) > 1E-6)
                        {
                            basis[k] = i;
                            newColumn = col;
                            replaced = true;
                            break;
                        }
                    }
                }
            }
        }
    }
}