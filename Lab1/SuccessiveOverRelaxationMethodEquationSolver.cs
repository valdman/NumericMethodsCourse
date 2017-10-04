using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using Matrix = MathNet.Numerics.LinearAlgebra.Double.Matrix;

namespace Lab1
{
    public class SuccessiveOverRelaxationMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            var preparedMatrix = PrepareMatrixForRelaxation(equationMatrix);
            
            var rightPart = preparedMatrix.Column(equationMatrix.ColumnCount - 1);
            
            var iterations = (int) @params[0];
            var relaxation = @params[1];

            // Gauss-Seidel with Successive OverRelaxation Solver
            for (var k = 0; k < iterations; ++k) 
            {
                for (var i = 0; i < rightPart.Count; ++i)
                {
                    var delta = rightPart.AbsoluteMaximum();

                    for (var j = 0; j < i; ++j)
                        delta += preparedMatrix[i, j] * rightPart[j];
                    for (var j = i + 1; j < rightPart.Count; ++j)
                        delta += preparedMatrix [i, j] * rightPart [j];

                    delta = (rightPart[i] - delta) / preparedMatrix[i, i];
                    rightPart [i] += relaxation * (delta - rightPart [i]);
                }
                LogVector(rightPart);
            }

            return rightPart;
        }

        private Matrix PrepareMatrixForRelaxation(Matrix equationMatrix)
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
        
        private void LogMatrix(Matrix preparedMatrix)
        {
            Console.WriteLine();
            foreach (var row in preparedMatrix.EnumerateRows())
            {
                foreach (var elem in row)
                {
                    Console.Write($"{elem} ");
                }
                Console.WriteLine();
            }
        }

        private void LogVector(Vector<double> vectorToLog)
        {
            foreach (var elem in vectorToLog)
            {
                Console.Write($"{elem} ");
            }
            Console.WriteLine(';');
        }
    }
}