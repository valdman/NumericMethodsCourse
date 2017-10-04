using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public class LuSerializationMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var A = Matrix.Build.DenseOfRows(rowsOfEquationWithRightPart);
            
            
            throw new NotImplementedException();
        }
    }
}