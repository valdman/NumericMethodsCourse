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
            var equationRows = new[]
            {
                new[] {1.0, 2, 3, 5},
                new[] {1.0, -1, 1, 0},
                new[] {1.0, -2, -3, 2}
            };
            
            var result = _equationSystemSolver.Solve(equationRows, 100, 1E-6);

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
                Console.Write($"x{i} = {string.Format("{0:0.00}", root)}; ");
                ++i;
            }
        }
    }
}