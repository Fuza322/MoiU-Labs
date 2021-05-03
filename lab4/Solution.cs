using System;
using System.Globalization;
using System.Linq;

namespace SimplexMethod
{
    public static class Program
    {
        // Enter the initial conditions
        static double[] InputDouble() => Console.ReadLine().Split(' ').Select(x => double.Parse(x, CultureInfo.InvariantCulture)).ToArray();

        // Enter the initial conditions
        static int[] InputInt() => Console.ReadLine().Split(' ').Select(x => int.Parse(x)).ToArray();

        // Enter the initial conditions
        static double[][] InputMatrix(int N, int M)
        {
            var matrix = DecMatrix(N, M);
            for (int i = 0; i < N; i++)
            {
                var input = InputDouble();
                for (int j = 0; j < M; j++)
                {
                    matrix[i][j] = input[j];
                }
            }
            return matrix;
        }

        // Geting a column of matrix
        static double[] Get_Column(double[][] matrix, int i)
        {
            var N = matrix.LengthRows();
            var result = new double[N];
            for (int j = 0; j < N; j++)
            {
                result[j] = matrix[j][i];
            }
            return result;
        }

        // Geting a row of matrix
        static double[] Get_Row(double[][] matrix, int i) => matrix[i];

        // Rows count calculation
        public static int LengthRows(this double[][] matrix)
        {
            return matrix.Length;
        }
        // Columns count calculation
        public static int LengthColumns(this double[][] matrix)
        {
            return matrix[0].Length;
        }

        // Multiply the vector by the matrix
        static double[] MulVectorOnMatrix(double[] vect, double[][] matrix)
        {
            var result = new double[matrix.LengthColumns()];
            for (int i = 0; i < matrix.LengthColumns(); i++)
            {
                result[i] = 0;
                for (int j = 0; j < matrix.LengthRows(); j++)
                {
                    result[i] += matrix[j][i] * vect[j];
                }
            }
            return result;
        }

        // Multiply the matrix by the vector
        static double[] MulMatrixOnVector(double[][] matrix, double[] vect)
        {
            var result = new double[matrix.LengthRows()];
            for (int i = 0; i < matrix.LengthRows(); i++)
            {
                result[i] = 0;
                for (int j = 0; j < matrix.LengthColumns(); j++)
                {
                    result[i] += matrix[i][j] * vect[j];
                }
            }
            return result;
        }

        // Vector multiplication
        static double MulVectors(double[] vect_A, double[] vect_B)
        {
            double result = 0;
            for (int i = 0; i < vect_A.Length; i++)
            {
                result += vect_A[i] * vect_B[i];
            }
            return result;
        }

        // Matrix Declaration
        public static double[][] DecMatrix(int N, int M)
        {
            var matrix = new double[N][];
            for (int i = 0; i < N; i++)
            {
                matrix[i] = new double[M];
            }
            return matrix;
        }

        // Sherman-Morrison's formul
        public static double[][] Sherman_morrison(double[][] matrix, double[] col, int colIndex)
        {
            var N = matrix.LengthRows();
            var vect1_l = MulMatrixOnVector(matrix, col);
            var vect2_l = new double[N];
            vect1_l.CopyTo(vect2_l, 0);
            vect2_l[colIndex] = -1;
            var matrix_Q = DecMatrix(N, N);
            for (int i = 0; i < N; i++)
            {
                matrix_Q[i][i] = 1;
                vect2_l[i] = -1 / vect1_l[colIndex] * vect2_l[i];
                matrix_Q[i][colIndex] = vect2_l[i];
            }
            //-------------------------------------------------------------------------------
            var matrix_C = DecMatrix(N, N);
            for (var i = 0; i < colIndex; i++)
            {
                for (var j = 0; j < N; j++)
                {
                    matrix_C[i][j] = matrix_Q[i][colIndex] * matrix[colIndex][j] + matrix_Q[i][i] * matrix[i][j];
                }
            }
            for (var j = 0; j < N; j++)
            {
                matrix_C[colIndex][j] = matrix_Q[colIndex][colIndex] * matrix[colIndex][j];
            }
            for (var i = colIndex + 1; i < N; i++)
            {
                for (var j = 0; j < N; j++)
                {
                    matrix_C[i][j] = matrix_Q[i][colIndex] * matrix[colIndex][j] + matrix_Q[i][i] * matrix[i][j];
                }
            }
            return matrix_C;
        }

        // Transpose matrix
        public static void InverseMatrix(double[][] matrix)
        {
            var N = matrix.LengthRows();
            for (int i = 1; i < N; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    var tmp = matrix[i][j];
                    matrix[i][j] = matrix[j][i];
                    matrix[j][i] = tmp;
                }
            }
        }

        // The formation of the basis matrix
        static double[][] CalculateBasisMatrix(double[][] matrix, int[] basis)
        {
            int N = basis.Length;
            var basis_Matrix = DecMatrix(N, N);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    basis_Matrix[i][j] = matrix[i][basis[j]];
                }
            }
            return basis_Matrix;
        }

        // Getting the basis vector
        static double[] CaltulateBasisVector(double[] vect, int[] basis) => basis.Select(e => vect[e]).ToArray();

        // Geting the Delta
        static double[] Get_Delta(double[] vect_u, double[][] matrix_A, double[] vect_c, int[] noNbasis)
        {
            var result = new double[matrix_A.LengthColumns()];
            foreach (var i in noNbasis)
            {
                result[i] = 0;
                for (int j = 0; j < matrix_A.LengthRows(); j++)
                {
                    result[i] += matrix_A[j][i] * vect_u[j];
                }
                result[i] -= vect_c[i];
            }
            return result;
        }

        // Getting the nobasis vector
        static int[] CalculateNonBasis(int[] basis, int N) => Enumerable.Range(0, N).Except(basis).ToArray();

        // Deltas rating
        static int? CheckDelta(double[] Δ, int[] noNbasis)
        {
            foreach (var i in noNbasis)
            {
                if (Δ[i] < 0)
                { return i; }
            }
            return null;
        }

        // Getting the Z-vector
        static double[] GetVect_Z(double[][] matrix_B, double[][] matrix_A, int badIndex)
        {
            var N = matrix_B.LengthRows();
            var result = new double[N];
            for (int i = 0; i < N; i++)
            {
                result[i] = 0;
                for (int j = 0; j < N; j++)
                {
                    result[i] += matrix_B[i][j] * matrix_A[j][badIndex];
                }
            }
            return result;
        }

        // Finding Teta 
        static double[] GetTeta(double[] vect_x, int[] basis, double[] vect_z)
        {
            var N = vect_z.Length;
            var teta = new double[N];
            for (int i = 0; i < N; i++)
            {
                teta[i] = vect_z[i] > 0 ? vect_x[basis[i]] / vect_z[i] : double.PositiveInfinity;
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
        static (double[], int[]) FindNewPlan(double[] oldPlan, int[] oldBasis, double min_teta, int minIndex, int J_0, double[] vect_z)
        {
            var N = oldPlan.Length;
            var M = oldBasis.Length;
            var newPlan = new double[N];
            newPlan[J_0] = min_teta;
            for (int i = 0; i < M; i++)
            {
                newPlan[oldBasis[i]] = oldPlan[oldBasis[i]] - min_teta * vect_z[i];
            }
            var newBasis = new int[M];
            oldBasis.CopyTo(newBasis, 0);
            newBasis[minIndex] = J_0;
            return (newPlan, newBasis);
        }

        // Geting the new matrix after removing the K-row and K-column from the initial solutions
        public static double[][] GetNewMatrix(double[][] matrix, int rowToDelete, int colToDelete)
        {
            var newMatrix = DecMatrix(matrix.LengthRows() - 1, matrix.LengthColumns() - 1);
            int rowIndex = 0;
            int colIndex = 0;
            for (int i = 0; rowIndex < matrix.LengthRows() - 1 && i < newMatrix.LengthRows(); i++)
            {
                rowIndex = i >= rowToDelete ? i + 1 : i;
                for (int j = 0; colIndex < matrix.LengthColumns() - 1 && j < newMatrix.LengthColumns(); j++)
                {
                    colIndex = j >= colToDelete ? j + 1 : j;
                    newMatrix[i][j] = matrix[rowIndex][colIndex];
                    colIndex = 0;
                }
            }
            return newMatrix;
        }

        // Algorithm of the main phase of the simplex-method
        static (double[], int[], double[][]) MainPhase(double[][] matrix_A, double[] vect_c, double[] vect_x, int[] basis, int N, double[][] matrix_B = null)
        {
            double[] vect_u = null;
            int? J_0 = 0;
            int min_J = 0;
            double min_teta = 0;
            if (matrix_B == null)
            {
                matrix_B = CalculateBasisMatrix(matrix_A, basis);
                vect_u = CaltulateBasisVector(vect_c, basis);
            }
            else
            {
                vect_u = MulVectorOnMatrix(CaltulateBasisVector(vect_c, basis), matrix_B);
            }
            while (true)
            {
                var noNbasis = CalculateNonBasis(basis, N);
                var Δ = Get_Delta(vect_u, matrix_A, vect_c, noNbasis);
                J_0 = CheckDelta(Δ, noNbasis);
                if (J_0 == null)
                {
                    return (vect_x, basis, matrix_B);
                }
                var vect_z = GetVect_Z(matrix_B, matrix_A, J_0.Value);
                if (vect_z.All(p => p < 1E-6))
                { return (null, null, null); }
                var teta = GetTeta(vect_x, basis, vect_z);
                (min_teta, min_J) = GetTeta_min(teta);
                (vect_x, basis) = FindNewPlan(vect_x, basis, min_teta, min_J, J_0.Value, vect_z);
                matrix_B = Sherman_morrison(matrix_B, Get_Column(matrix_A, J_0.Value), min_J);
                vect_u = MulVectorOnMatrix(CaltulateBasisVector(vect_c, basis), matrix_B);
            }
        }

        // Algorthm of the first phase simplex-method
        public static (double[], int[], double[][], double[][]) FirstPhase(double[][] matrix_A, double[] vect_b, double[] vect_c, int N, int M)
        {
            var matrixExpanded = DecMatrix(M, N + M);
            vect_c = new double[N + M];
            for (int i = 0; i < M; i++)
            {
                if (vect_b[i] < 0)
                {
                    vect_b[i] = -vect_b[i];
                    for (int j = 0; j < N; j++)
                    {
                        matrix_A[i][j] = -matrix_A[i][j];
                    }
                }
                for (int j = 0; j < N; j++)
                {
                    matrixExpanded[i][j] = matrix_A[i][j];
                }
            }
            var vect_x = new double[N + M];
            var basis = new int[M];
            for (int i = 0; i < M; i++)
            {
                vect_c[i + N] = -1;
                matrixExpanded[i][i + N] = 1;
                vect_x[i + N] = vect_b[i];
                basis[i] = i + N;
            }
            double[][] rev_Ab = null;
            (vect_x, basis, rev_Ab) = MainPhase(matrixExpanded, vect_c, vect_x, basis, N + M);
            if (vect_x == null)
            {
                Console.WriteLine("Unbounded");
                return (null, null, null, null);
            }
            //
            for (int i = 0; i < M; i++)
            {
                if (Math.Abs(vect_x[i + N]) > 1E-6)
                {
                    Console.WriteLine("No solution");
                    return (null, null, null, null);
                }
            }
            var replaced = false;
            int k = 0;
            while (true)
            {
                // 
                if (basis.All(e => e < N))
                {
                    vect_x = vect_x.Take(N).ToArray();

                    var newMatrix_A = DecMatrix(M, N);

                    for (int i = 0; i < M; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            newMatrix_A[i][j] = matrixExpanded[i][j];
                        }
                    }
                    return (vect_x, basis, newMatrix_A, rev_Ab);
                }
                // 
                var Ab = CalculateBasisMatrix(matrixExpanded, basis);
                var noNbasis = CalculateNonBasis(basis, N);
                k = 0;
                for (int i = 0; i < M; i++)
                {
                    if (basis[i] > N)
                    {
                        k = i;
                        break;
                    }
                }
                replaced = false;
                var lessThenM = noNbasis.Where(i => i < M);
                foreach (var i in lessThenM)
                {
                    var col = Get_Column(matrixExpanded, i);
                    var alpha = MulVectors(Get_Row(rev_Ab, k), col);
                    if (Math.Abs(alpha) > 1E-6)
                    {
                        basis[k] = i;
                        replaced = true;
                        rev_Ab = Sherman_morrison(rev_Ab, col, k);
                        break;
                    }
                }
                if (!replaced)
                {
                    matrixExpanded = GetNewMatrix(matrixExpanded, basis[k] - N, basis[k]);
                    var newBasis = basis.Select(i => i > basis[k] ? i - 1 : i).ToList();
                    M--;
                    newBasis.RemoveAt(k);
                    basis = newBasis.ToArray();

                    // Finding a new inversion after we throw out the row and column
                    var e = new double[rev_Ab.LengthRows()];
                    e[k] = 1;
                    var matrix = Sherman_morrison(rev_Ab, e, k);
                    InverseMatrix(matrix);
                    matrix = Sherman_morrison(rev_Ab, e, k);
                    var newMatrix = GetNewMatrix(matrix, k, k);
                    rev_Ab = newMatrix;
                }
            }
        }

        static void Main()
        {
            var input = InputInt();
            var M = input[0];
            var N = input[1];
            var matrix_A = InputMatrix(M, N);
            var vect_b = InputDouble();
            var vect_c = InputDouble();

            (var x, var basis, var newA, var B) = FirstPhase(matrix_A, vect_b, vect_c, N, M);
            if (x == null)
            {
                return;
            }
            var result = MainPhase(newA, vect_c, x, basis, N, B).Item1;
            if (result == null)
            {
                Console.WriteLine("Unbounded");
                return;
            }
            Console.WriteLine("Bounded");
            // Result output
            for (var i = 0; i < result.Length; i++)
            {
                Console.Write(result[i].ToString(CultureInfo.GetCultureInfo("en-US")) + " ");
            }
        }
    }
}