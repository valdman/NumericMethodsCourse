using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public static class MatrixExtensions
    {
        public static Matrix PrepareMatrixForIterations(this Matrix<double> equationMatrix)
        {
            var preparedMatrix = equationMatrix.Clone();
            
            for (var i = 0; i < equationMatrix.RowCount; i++)
            {
                for (var j = 0; j < equationMatrix.ColumnCount - 1; j++)
                {
                    preparedMatrix[i, j] /=  -equationMatrix[i, i];
                }

                preparedMatrix[i, equationMatrix.ColumnCount - 1] /= equationMatrix[i, i];
            }

            //LogMatrix();
            return preparedMatrix as Matrix;
        }

        public static (Matrix L, Matrix U) LuSerialize(this Matrix<double> equationMatrix)
        {
            var n = equationMatrix.RowCount;
            var m = equationMatrix.ColumnCount;
            if (n != m)
                throw new ArgumentException("For LU method coefficient matrix must be square");

            var l = new DenseMatrix(n);
            var u = new DenseMatrix(n);
            for (var i = 0; i < n; i++)
            {
                for (var j = 0; j < n; j++)
                {
                    u[0, i] = equationMatrix[0, i];
                    l[i, 0] = equationMatrix[i, 0] / u[0, 0];
                    double sum = 0;
                    for (var k = 0; k < i; k++)
                    {
                        sum += l[i, k] * u[k, j];
                    }
                    u[i, j] = equationMatrix[i, j] - sum;
                    if (i > j)
                    {
                        l[j, i] = 0;
                    }
                    else
                    {
                        sum = 0;
                        for (var k = 0; k < i; k++)
                        {
                            sum += l[j, k] * u[k, i];
                        }
                        l[j, i] = (equationMatrix[j, i] - sum) / u[i, i];
                    }
                }
            }
            
            return (l, u);
        }
    }
}