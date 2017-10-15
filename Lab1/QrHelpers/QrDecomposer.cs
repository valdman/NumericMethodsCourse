using System;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace Lab1.QrHelpers
{
    public class QrDecomposition
    {
        public static UserQR Create(Matrix<double> matrix, QRMethod method = QRMethod.Full)
        {
            if (matrix.RowCount < matrix.ColumnCount)
                throw new ArgumentException("Dimensions doesn't match");
            var length = Math.Min(matrix.RowCount, matrix.ColumnCount);
            var numArray = new double[length][];
            Matrix<double> matrix1;
            Matrix<double> matrix2;
            if (method == QRMethod.Full)
            {
                matrix1 = matrix.Clone();
                matrix2 = Matrix<double>.Build.SameAs(matrix, matrix.RowCount, matrix.RowCount, true);
                for (var index = 0; index < matrix.RowCount; ++index)
                    matrix2.At(index, index, 1.0);
                for (var index = 0; index < length; ++index)
                {
                    numArray[index] = UserQR.GenerateColumn(matrix1, index, index);
                    UserQR.ComputeQR(numArray[index], matrix1, index, matrix.RowCount, index + 1, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
                }
                for (var index = length - 1; index >= 0; --index)
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index, matrix.RowCount,
                        Control.MaxDegreeOfParallelism);
            }
            else
            {
                matrix2 = matrix.Clone();
                for (var index = 0; index < length; ++index)
                {
                    numArray[index] = UserQR.GenerateColumn(matrix2, index, index);
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index + 1, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
                }
                matrix1 = matrix2.SubMatrix(0, matrix.ColumnCount, 0, matrix.ColumnCount);
                matrix2.Clear();
                for (var index = 0; index < matrix.ColumnCount; ++index)
                    matrix2.At(index, index, 1.0);
                for (var index = length - 1; index >= 0; --index)
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
            }
            return new UserQR(matrix2, matrix1, method);
        }
    }
}