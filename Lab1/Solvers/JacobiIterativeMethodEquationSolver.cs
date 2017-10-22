using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1.Solvers
{
    public class JacobIIterativeMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);

            var n = equationMatrix.ColumnCount - 1;
            var rightPart = equationMatrix.Column(n);
            var leftMatrix = equationMatrix.RemoveColumn(n);

            if(leftMatrix.RowCount != leftMatrix.ColumnCount) throw new ArgumentException("For Simple Jacobi method matrix must be square");

            //if (!TestConvergenceForJacobi(leftMatrix)) throw new ArgumentException("Method does not converge on this matrix");

            var (d, l, u) = leftMatrix.DluDecompose();
            var b = d.Inverse() * (d - leftMatrix);
            var g = d.Inverse() * rightPart;

            var x = rightPart.Clone();

            var iterations = (int) @params[0];

            for (var k = 0; k < iterations; k++)
            {
                var y = b * x + g;
                x = y.Clone();
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