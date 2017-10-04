using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public class SuccessiveOverRelaxationMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            var preparedMatrix = equationMatrix.Clone();
            
            for (var i = 0; i < equationMatrix.RowCount; i++)
            {
                for (var j = 0; j < equationMatrix.ColumnCount; j++)
                {
                    preparedMatrix[i, j] /=  -equationMatrix[i, i];
                }
            }
            
            var rightPart = equationMatrix.Column(equationMatrix.ColumnCount - 1);
            
            var iterations = (int) @params[0];
            var relaxation = @params[1];

            // Gauss-Seidel with Successive OverRelaxation Solver
            for (var k = 0; k < iterations; ++k) 
            {
                for (var i = 0; i < rightPart.Count; ++i) 
                {
                    double delta = 0.0f;

                    for (var j = 0; j < i; ++j)
                        delta += equationMatrix[ i, j] * rightPart[j];
                    for (var j = i + 1; j < rightPart.Count; ++j)
                        delta += equationMatrix [i, j] * rightPart [j];

                    delta = (rightPart[i] - delta) / equationMatrix[i, i];
                    rightPart [i] += relaxation * (delta - rightPart [i]);
                }
            }

            return rightPart;
        }
    }
}