using System.Web.Http;
using Newtonsoft.Json;
using Newtonsoft.Json.Serialization;
using Owin;

namespace DWSIM.Api
{
    public class Startup
    {
        public void Configuration(IAppBuilder app)
        {
            var config = new HttpConfiguration();
            config.MapHttpAttributeRoutes();

            // camelCase JSON, ignore reference loops (DWSIM object graphs are cyclic).
            var json = config.Formatters.JsonFormatter.SerializerSettings;
            json.ContractResolver = new CamelCasePropertyNamesContractResolver();
            json.ReferenceLoopHandling = ReferenceLoopHandling.Ignore;

            app.UseWebApi(config);
        }
    }
}
