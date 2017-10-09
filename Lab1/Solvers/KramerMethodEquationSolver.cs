using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1.Solvers
{
    public class KramerMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquation, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquation);

            var coefficientMatrix = equationMatrix.SubMatrix(0, equationMatrix.RowCount, 0, equationMatrix.ColumnCount - 1);
            var mainDeterminant = coefficientMatrix.Determinant();

            var rightPartVector = equationMatrix.Column(equationMatrix.ColumnCount - 1);

            var answerVector = new double[coefficientMatrix.ColumnCount];

            for (var i = 0; i < coefficientMatrix.ColumnCount; i++)
            {
                var currentMatrix = coefficientMatrix.Clone();
                currentMatrix.SetColumn(i, rightPartVector);

                answerVector[i] = currentMatrix.Determinant() / mainDeterminant;
            }

            return answerVector;
        }
    }
}