using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public class QrDecompositionMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            var solver = new QrDecomposer(equationMatrix);

            var (q, r) = (solver.OrthogonalFactor, solver.UpperTriangularFactor);

            if (Equals(q * r, equationMatrix))
            {
                Console.WriteLine("QR decomposing works");
            }
            
            throw new NotImplementedException();
        }
    }
}