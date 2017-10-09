using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Properties;

namespace Lab1.QrHelpers
{
    public sealed class UserQR : QR
    {
        public static UserQR Create(Matrix<double> matrix, QRMethod method = QRMethod.Full)
        {
            if (matrix.RowCount < matrix.ColumnCount)
                throw new ArgumentException("Dimensitons doesn't match");
            int length = Math.Min(matrix.RowCount, matrix.ColumnCount);
            double[][] numArray = new double[length][];
            Matrix<double> matrix1;
            Matrix<double> matrix2;
            if (method == QRMethod.Full)
            {
                matrix1 = matrix.Clone();
                matrix2 = Matrix<double>.Build.SameAs<double>(matrix, matrix.RowCount, matrix.RowCount, true);
                for (int index = 0; index < matrix.RowCount; ++index)
                    matrix2.At(index, index, 1.0);
                for (int index = 0; index < length; ++index)
                {
                    numArray[index] = UserQR.GenerateColumn(matrix1, index, index);
                    UserQR.ComputeQR(numArray[index], matrix1, index, matrix.RowCount, index + 1, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
                }
                for (int index = length - 1; index >= 0; --index)
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index, matrix.RowCount,
                        Control.MaxDegreeOfParallelism);
            }
            else
            {
                matrix2 = matrix.Clone();
                for (int index = 0; index < length; ++index)
                {
                    numArray[index] = UserQR.GenerateColumn(matrix2, index, index);
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index + 1, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
                }
                matrix1 = matrix2.SubMatrix(0, matrix.ColumnCount, 0, matrix.ColumnCount);
                matrix2.Clear();
                for (int index = 0; index < matrix.ColumnCount; ++index)
                    matrix2.At(index, index, 1.0);
                for (int index = length - 1; index >= 0; --index)
                    UserQR.ComputeQR(numArray[index], matrix2, index, matrix.RowCount, index, matrix.ColumnCount,
                        Control.MaxDegreeOfParallelism);
            }
            return new UserQR(matrix2, matrix1, method);
        }

        internal UserQR(Matrix<double> q, Matrix<double> rFull, QRMethod method)
            : base(q, rFull, method)
        {
        }

        internal static double[] GenerateColumn(Matrix<double> a, int row, int column)
        {
            int length = a.RowCount - row;
            double[] numArray = new double[length];
            for (int row1 = row; row1 < a.RowCount; ++row1)
            {
                numArray[row1 - row] = a.At(row1, row);
                a.At(row1, row, 0.0);
            }
            double num1 = Math.Sqrt(((IEnumerable<double>) numArray).Sum<double>((Func<double, double>) (t => t * t)));
            if (row == a.RowCount - 1 || num1 == 0.0)
            {
                a.At(row, column, -numArray[0]);
                numArray[0] = 1.4142135623731;
                return numArray;
            }
            double num2 = 1.0 / num1;
            if (numArray[0] < 0.0)
                num2 *= -1.0;
            a.At(row, column, -1.0 / num2);
            for (int index = 0; index < length; ++index)
                numArray[index] *= num2;
            ++numArray[0];
            double num3 = Math.Sqrt(1.0 / numArray[0]);
            for (int index = 0; index < length; ++index)
                numArray[index] *= num3;
            return numArray;
        }

        internal static void ComputeQR(double[] u, Matrix<double> a, int rowStart, int rowDim, int columnStart,
            int columnDim, int availableCores)
        {
            if (rowDim < rowStart || columnDim < columnStart)
                return;
            int num1 = columnDim - columnStart;
            for (int column = columnStart; column < columnDim; ++column)
            {
                double num2 = 0.0;
                for (int row = rowStart; row < rowDim; ++row)
                    num2 += u[row - rowStart] * a.At(row, column);
                for (int row = rowStart; row < rowDim; ++row)
                    a.At(row, column, a.At(row, column) - u[row - rowStart] * num2);
            }
        }

        public override void Solve(Matrix<double> input, Matrix<double> result)
        {
            if (input.ColumnCount != result.ColumnCount)
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            if (this.FullR.RowCount != input.RowCount)
                throw new ArgumentException(Resources.ArgumentMatrixSameRowDimension);
            if (this.FullR.ColumnCount != result.RowCount)
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            Matrix<double> matrix = input.Clone();
            double[] numArray = new double[this.FullR.RowCount];
            for (int column = 0; column < input.ColumnCount; ++column)
            {
                for (int row = 0; row < this.FullR.RowCount; ++row)
                    numArray[row] = matrix.At(row, column);
                for (int index = 0; index < this.FullR.RowCount; ++index)
                {
                    double num = 0.0;
                    for (int row = 0; row < this.FullR.RowCount; ++row)
                        num += this.Q.At(row, index) * numArray[row];
                    matrix.At(index, column, num);
                }
            }
            for (int index = this.FullR.ColumnCount - 1; index >= 0; --index)
            {
                for (int column = 0; column < input.ColumnCount; ++column)
                    matrix.At(index, column, matrix.At(index, column) / this.FullR.At(index, index));
                for (int row = 0; row < index; ++row)
                {
                    for (int column = 0; column < input.ColumnCount; ++column)
                        matrix.At(row, column,
                            matrix.At(row, column) - matrix.At(index, column) * this.FullR.At(row, index));
                }
            }
            for (int row = 0; row < this.FullR.ColumnCount; ++row)
            {
                for (int column = 0; column < matrix.ColumnCount; ++column)
                    result.At(row, column, matrix.At(row, column));
            }
        }

        public override void Solve(Vector<double> input, Vector<double> result)
        {
            if (this.FullR.RowCount != input.Count)
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            if (this.FullR.ColumnCount != result.Count)
                throw new ArgumentException("Dimensitons doesn't match");
            Vector<double> vector1 = input.Clone();
            double[] numArray = new double[this.FullR.RowCount];
            for (int index = 0; index < this.FullR.RowCount; ++index)
                numArray[index] = vector1[index];
            for (int column = 0; column < this.FullR.RowCount; ++column)
            {
                double num = 0.0;
                for (int row = 0; row < this.FullR.RowCount; ++row)
                    num += this.Q.At(row, column) * numArray[row];
                vector1[column] = num;
            }
            for (int index1 = this.FullR.ColumnCount - 1; index1 >= 0; --index1)
            {
                Vector<double> vector2 = vector1;
                int index2 = index1;
                vector2[index2] = vector2[index2] / this.FullR.At(index1, index1);
                for (int row = 0; row < index1; ++row)
                {
                    Vector<double> vector3 = vector1;
                    int index3 = row;
                    vector3[index3] = vector3[index3] - vector1[index1] * this.FullR.At(row, index1);
                }
            }
            for (int index = 0; index < this.FullR.ColumnCount; ++index)
                result[index] = vector1[index];
        }
    }
}