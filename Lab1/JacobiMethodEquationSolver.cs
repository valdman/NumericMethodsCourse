using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public class JacobiMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            var n = equationMatrix.ColumnCount - 1;
            
            var rightPart = equationMatrix.Column(n);
            var a = equationMatrix.RemoveColumn(n).PrepareMatrixForIterations();

            var iterations = (int) @params[0];

            var x = new double[n];
            for (var q = 0; q < iterations; q++)
            {
                for (var i = 0; i < rightPart.Count; i++)
                {
                    var sum = 0.0;
                    for (var j = 0; j < n; j++)
                    {
                        if(j == i) continue;
                        sum += (a[i, j] / a[i, i]) * x[j];
                    }

                    x[i] = rightPart[i] / a[i, i] - sum;
                }
            }

            return x;
        }
    }
}