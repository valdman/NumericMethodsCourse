using System;
using System.Collections.Generic;
using System.Linq;
using Lab1;
using Microsoft.Extensions.DependencyInjection;

namespace NumercMethods
{
    public class Worker
    {
        private readonly IEquationSystemSolver _equationSystemSolver;

        public Worker(IEquationSystemSolver equationSystemSolver)
        {
            _equationSystemSolver = equationSystemSolver;
        }

        public void Job()
        {
            //var n = ReadNumberOfRows();
            //var equationRows = ReadEquationRows(n);
            var oldEquationRows = new[]
            {
                new[] {7.0, 1, 2, 2, 1},
                new[] {0.0, 5, 1, 1, -3},
                new[] {2.0, 1, 8, 1, 0},
                new[] {2.0, 2, 1, 6, 4}
            };

            var equationRows = new[]
            {
                new double[] {2, 1, 1},
                new double[] {1, 2, -1}
            };
            
            var result = _equationSystemSolver.Solve(oldEquationRows, 1000, 1E-5);

            WriteVectorAsRoots(result);
        }

        private int ReadNumberOfRows()
        {
            Console.WriteLine("Write number of rows");
            return int.Parse(Console.ReadLine());
        }

        private IEnumerable<double[]> ReadEquationRows(int numberOfRows)
        {
            Console.WriteLine("Writel equation rows coefficients (with right part)");
            var equationRows = new List<double[]>();
            
            for (var i = 0; i < numberOfRows; i++)
            {
                var coefficientsForCurrentRow = Console.ReadLine()
                    .Split(' ', StringSplitOptions.RemoveEmptyEntries)
                    .Select(double.Parse).ToArray();
                equationRows.Add(coefficientsForCurrentRow);
            }

            return equationRows;
        }

        private void WriteVectorAsRoots(IEnumerable<double> roots)
        {
            Console.WriteLine();
            var i = 0;
            foreach (var root in roots)
            {
                Console.Write($"x{i} = {string.Format("{0:0.0000}", root)}; ");
                ++i;
            }
        }
    }
}