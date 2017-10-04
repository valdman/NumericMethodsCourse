using System.Collections.Generic;

namespace Lab1
{
    public interface IEquationSystemSolver
    {
        IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params);
    }
}