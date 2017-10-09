using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1.Solvers
{
    public class JacobiIterativeMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);

            var n = equationMatrix.ColumnCount - 1;
            var rightPart = equationMatrix.Column(n);
            var leftMatrix = equationMatrix.RemoveColumn(n);

            if(leftMatrix.RowCount != leftMatrix.ColumnCount) throw new ArgumentException("For Simple Jacobi method matrix must be square");

            if (!TestConvergenceForJacobi(leftMatrix)) throw new ArgumentException("Method does not converge on this matrix");

            var x = new double[n];
            var y = new double[n];

            var iterations = (int) @params[0];

            for (var k = 0; k < iterations; k++)
            {
                for (var i = 0; i < n; i++)
                {
                    var sum = 0.0;
                    for (var j = 0; j < n; j++)
                    {
                        if(j == i) continue;
                        sum += leftMatrix[i, j] * x[j];
                    }
                    /*
                     * Using this we will get Gauss-Zeidel method
                    x[i] = (rightPart[i] - sum) / leftMatrix[i, i];
                    */
                    y[i] = (rightPart[i] - sum) / leftMatrix[i, i];
                }
                x = y;
            }

            return x;
        }

        private bool TestConvergenceForJacobi(Matrix<double> leftMatrix)
        {
            var (d, l, u) = leftMatrix.DluDecompose();
            var convergentor = d.Inverse() * (l + u);
            var convergentorNorm = convergentor.InfinityNorm();
            return (convergentorNorm < 1);
        }
    }
}