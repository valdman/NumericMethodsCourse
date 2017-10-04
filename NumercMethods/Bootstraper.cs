using Lab1;
using Microsoft.Extensions.DependencyInjection;

namespace NumercMethods
{
    public class Bootstraper
    {
        private readonly ServiceCollection _service;

        public Bootstraper(ServiceCollection service)
        {
            _service = service;
        }

        public void Configure()
        {
            _service
                .AddLogging()
                .AddSingleton<IEquationSystemSolver, KramerMethodEquationSolver>()
                .BuildServiceProvider();
        }
    }
}