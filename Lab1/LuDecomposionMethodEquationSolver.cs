﻿using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Lab1
{
    public class LuDecomposionMethodEquationSolver : IEquationSystemSolver
    {
        public IEnumerable<double> Solve(IEnumerable<double[]> rowsOfEquationWithRightPart, params double[] @params)
        {
            var equationMatrix = DenseMatrix.OfRows(rowsOfEquationWithRightPart);
            var n = equationMatrix.ColumnCount - 1;
            
            var rightPart = equationMatrix.Column(n);
            var leftMatrix = equationMatrix.RemoveColumn(n);

            var (l, u) = leftMatrix.LuSerialize();
            var y = new double[n];
            var x = new double[n];

            for (int i = 0; i < n; i++)
            {
                var sumForX = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sumForX += l[i, j] * y[j];
                }
                y[i] = (rightPart[i] - sumForX) / l[i, i];
            }

            for (int i = n - 1; i >= 0; i--)
            {
                var sumForY = 0.0;
                for (int j = n - 1; j > i; j--)
                {
                    sumForY += u[i, j] * x[j];
                }
                x[i] = (y[i] - sumForY) / u[i, i];
            }

            return x;
        }
    }
}