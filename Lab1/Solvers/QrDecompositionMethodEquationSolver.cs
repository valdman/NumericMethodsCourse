using System;
using System.Collections.Generic;
using Lab1.QrHelpers;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1.Solvers
{
    public class QrDecompositionMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            
            var n = equationMatrix.ColumnCount - 1;
            var rightPart = equationMatrix.Column(n);
            var leftPart = equationMatrix.RemoveColumn(n);
                
            var qr = QrDecomposition.Create(leftPart);
            LogMatrix(qr.Q);
            LogMatrix(qr.R);
            Console.WriteLine();

            return qr.Solve(rightPart);
        }
        
        private void LogMatrix(Matrix<double> preparedMatrix)
        {
            Console.WriteLine();
            foreach (var row in preparedMatrix.EnumerateRows())
            {
                foreach (var elem in row)
                {
                    Console.Write($"{string.Format("{0:0.0000}", elem)} ");
                }
                Console.WriteLine();
            }
        }
    }
}