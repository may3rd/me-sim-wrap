using System;
using Microsoft.Owin.Hosting;

namespace DWSIM.Api
{
    // Self-hosted Web API console. Mirrors DWSIM.Apps.TCPServer's console-exe model — no IIS.
    class Program
    {
        static void Main(string[] args)
        {
            var baseUrl = Environment.GetEnvironmentVariable("DWSIM_API_URL") ?? "http://localhost:9000";

            // Build the engine once at startup: property packages etc. load here (slow — do it before serving).
            Console.WriteLine("Initializing DWSIM engine...");
            FlowsheetStore.Init();
            Console.WriteLine("Engine ready.");

            using (WebApp.Start<Startup>(baseUrl))
            {
                Console.WriteLine("DWSIM API listening on " + baseUrl);
                Console.WriteLine("Press Enter to stop.");
                Console.ReadLine();
            }
        }
    }
}
