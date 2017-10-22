using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using Matrix = MathNet.Numerics.LinearAlgebra.Double.Matrix;

namespace Lab1.Solvers
{
    public class SuccessiveOverRelaxationMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {      
            var relaxationCoefficient = @params[0];
            var eps = @params[1];
            
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            
            var b = equationMatrix.Column(equationMatrix.ColumnCount - 1);
            var a = equationMatrix.RemoveColumn(equationMatrix.ColumnCount - 1);
            
            var n = a.RowCount;
            
            var cur = b.Clone();

            int counter = 0;
            
            while (true)
            {
                counter++;
                var old = cur.Clone();
                for (int i = 0;i < n; ++i) {
                    double sum = 0;
                    for (int j = 0; j < n; ++j) 
                    {
                        if (j == i) continue;
                        sum -= (a[i, j] * cur[j]);
                    }
                    sum += b[i];
                    sum *= (relaxationCoefficient / a[i, i]);
                    cur[i] = sum+(1.0 - relaxationCoefficient)*cur[i];
                }
                double mx = Math.Abs(cur[0] - old[0]);
                for (int i = 0; i < n; ++i) {
                    mx = Math.Max(Math.Abs(cur[i] - old[i]), mx);
                }
                if (mx < eps) break;
            }
            Console.WriteLine(counter);
            return cur;
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