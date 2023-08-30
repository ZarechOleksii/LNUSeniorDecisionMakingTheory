using MathNet.Numerics.LinearAlgebra;
using System;

namespace Lab5
{
    public static class Utils
    {
        public static double[,] Normalization(double[,] matrix)
        {
            int size1 = matrix.GetLength(0);
            int size2 = matrix.GetLength(1);

            for (int i = 0; i < size2; ++i)
            {
                double min = matrix[0, i];
                double max = matrix[0, i];
                for (int j = 0; j < size1; ++j)
                {
                    if (matrix[j, i] > max)
                        max = matrix[j, i];
                    if (matrix[j, i] < min)
                        min = matrix[j, i];
                }
                for (int j = 0; j < size1; ++j)
                {
                    matrix[j, i] = (matrix[j, i] - min) / (max - min);
                }
            }

            return matrix;
        }

        public static Matrix<double> SmallerSquareMethod(Matrix<double> matrix, Matrix<double> vector)
        {
            return (matrix.Transpose() * matrix).Inverse() * (matrix.Transpose() * vector);
        }

        #region GradientMethod
        public static Matrix<double> GradientMethod(Matrix<double> matrix, Matrix<double> vector)
        {
            var W = matrix.Transpose();
            var x0 = Matrix<double>.Build.Dense(matrix.RowCount, 1);
            for (int i = 0; i < matrix.RowCount; i++)
            {
                x0[i, 0] = vector[i, 0] / matrix[i, 0]; //matrix[i, i]
            }

            Matrix<double> x = null;

            do
            {
                x?.CopyTo(x0);

                var Rp = Mul1(matrix, x0);

                var R = Matrix<double>.Build.Dense(Rp.RowCount, 1);
                for (int i = 0; i < Rp.RowCount; i++)
                {
                    R[i, 0] = Rp[i, 0] - vector[i, 0];
                }

                var D3 = matrix * W;
                var D = Mul1(D3, R);
                var D1 = ScalarProduction(R, D);
                var D2 = ScalarProduction(D, D);
                var temp = D1 / D2;
                x = SubGradientMethod(x0, R, W, temp);
            }
            while (!Difference(x, x0));

            return x;
        }
        private static bool Difference(Matrix<double> Xn, Matrix<double> X0)
        {
            var res = false;
            for (int i = 0; i < Xn.RowCount; i++)
            {
                if (Math.Abs(Xn[i, 0] - X0[i, 0]) > 0.001)
                {
                    return true;
                }
            }
            return res;
        }
        private static Matrix<double> Mul1(Matrix<double> A, Matrix<double> X0)
        {
            var xr = Matrix<double>.Build.Dense(A.RowCount, 1);
            for (int i = 0; i < A.RowCount; i++)
            {
                var sum = 0.0;
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    sum += A[i, j] * X0[j, 0];
                }
                xr[i, 0] = sum;
            }

            return xr;
        }
        private static double ScalarProduction(Matrix<double> X, Matrix<double> Y)
        {
            var s = 0.0;
            for (var i = 0; i < X.RowCount; i++)
            {
                s += X[i, 0] * Y[i, 0];
            }
            return s;
        }
        private static Matrix<double> SubGradientMethod(Matrix<double> X0, Matrix<double> R, Matrix<double> W, double temp)
        {
            var V = Mul1(W, R);
            var xr = Matrix<double>.Build.Dense(X0.RowCount, 1);
            for (int i = 0; i < X0.RowCount; i++)
            {
                xr[i, 0] = X0[i, 0] - (temp * V[i, 0]);
            }
            return xr;
        }
        #endregion

        private static double GetT(double x, int m)
        {
            if (m == 0)
                return 1;
            if (m == 1)
                return x;
            return 2 * x * GetT(x, m - 1) - GetT(x, m - 2);
        }
        public static double GetTStar(double x, int m)
        {
            if (m == 0)
                return 0.5;
            return GetT(2 * x - 1, m);
        }

        public static double[] GetB(double[,] y)
        {
            int size1 = y.GetLength(0);
            int size2 = y.GetLength(1);

            double[] res = new double[size1];
            for (int i = 0; i < size1; ++i)
            {
                double min = y[i, 0];
                double max = y[i, 0];
                for (int j = 0; j < size2; ++j)
                {
                    if (y[i, j] > max)
                        max = y[i, j];
                    if (y[i, j] < min)
                        min = y[i, j];
                }
                res[i] = (max + min) / 2;
            }
            return res;
        }

        public static double[,] GetPsi(int size, double[,] Lambda, double[,] x)
        {
            int m = x.GetLength(1);
            double[,] res = new double[size, m];
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    double sum = 0;
                    for (int p = 0; p < Lambda.GetLength(1); ++p)
                    {
                        sum += Lambda[j, p] * Math.Log(1e-9 + 1 + GetTStar(x[i, j], p));
                    }
                    res[i, j] = Math.Exp(sum) - 1;
                }
            }
            return res;
        }

        public static double[,] GetA(double[,] psi, double[,] y)
        {
            Matrix<double> _psi = Matrix<double>.Build.DenseOfArray(psi);
            Matrix<double> _y = Matrix<double>.Build.DenseOfArray(y);

            for (int i = 0; i < y.GetLength(0); ++i)
            {
                for (int j = 0; j < y.GetLength(1); ++j)
                {
                    _y[i, j] = Math.Log(1 + 1e-9 + _y[i, j]);
                }
            }
            for (int i = 0; i < psi.GetLength(0); ++i)
            {
                for (int j = 0; j < psi.GetLength(1); ++j)
                {
                    _psi[i, j] = Math.Log(1e-9 + 1 + _psi[i, j]);
                }
            }

            Matrix<double> resA = SmallerSquareMethod(_psi, _y);

            double[,] res = new double[y.GetLength(1), psi.GetLength(1)];
            for (int i = 0; i < y.GetLength(1); ++i)
            {
                for (int j = 0; j < psi.GetLength(1); ++j)
                {
                    res[i, j] = resA[j, i];
                }
            }
            return res;
        }

        public static double[] GetFi(int size, double[,] A, double[,] Psi, int index)
        {
            double[] res = new double[size];
            for (int i = 0; i < size; ++i)
            {
                double sum = 0;
                for (int j = 0; j < A.GetLength(1); ++j)
                {
                    sum += A[index, j] * Math.Log(Psi[i, j] + 1);
                }
                res[i] = Math.Exp(sum) - 1; ;
            }
            return res;
        }

        public static double[] GetC(
            int size, double[,] Y, int index,
            double[,] A1, double[,] A2, double[,] A3,
            double[,] Psi1, double[,] Psi2, double[,] Psi3)
        {
            Matrix<double> _Y = Matrix<double>.Build.Dense(size, 1);
            for (int i = 0; i < size; ++i)
            {
                _Y[i, 0] = Math.Log(Y[i, index] + 1);
            }
            Matrix<double> matrix = Matrix<double>.Build.Dense(size, 3);
            for (int i = 0; i < size; ++i)
            {
                matrix[i, 0] = Math.Log(1e-9 + 1 + GetFi(size, A1, Psi1, index)[i]);
                matrix[i, 1] = Math.Log(1e-9 + 1 + GetFi(size, A2, Psi2, index)[i]);
                matrix[i, 2] = Math.Log(1e-9 + 1 + GetFi(size, A3, Psi3, index)[i]);
            }

            Matrix<double> resC = SmallerSquareMethod(matrix, _Y);

            double[] res = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                res[i] = resC[i, 0];
            }
            return res;
        }
    }
}
