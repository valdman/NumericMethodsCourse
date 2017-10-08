using Lab1;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Logging;

namespace NumercMethods
{
    class Program
    {
        static void Main(string[] args)
        {
            //setup our DI
            var serviceProvider = new ServiceCollection()
                .AddLogging()
                .AddSingleton<Worker>()
                .AddSingleton<IEquationSystemSolver, QrDecompositionMethodEquationSolver>()
                .BuildServiceProvider();

            serviceProvider
                .GetService<ILoggerFactory>();

            var logger = serviceProvider.GetService<ILoggerFactory>()
                .CreateLogger<Program>();
            logger.LogDebug("Starting application");

            //do the actual work here
            var worker = serviceProvider.GetService<Worker>();
            worker.Job();

            logger.LogDebug("All done!");
        }
    }
}