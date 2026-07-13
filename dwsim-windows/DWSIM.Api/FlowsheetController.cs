using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Threading.Tasks;
using System.Web.Http;
using DWSIM.Interfaces;
using DWSIM.Interfaces.Enums;

namespace DWSIM.Api
{
    [RoutePrefix("api/flowsheet")]
    public class FlowsheetController : ApiController
    {
        // POST /api/flowsheet  -> new empty flowsheet
        [HttpPost, Route("")]
        public IHttpActionResult CreateEmpty()
        {
            var id = FlowsheetStore.Create(auto => auto.CreateFlowsheet());
            return Ok(new { id });
        }

        // POST /api/flowsheet/upload  -> body is a raw .dwxmz file; returns the session id
        [HttpPost, Route("upload")]
        public async Task<IHttpActionResult> Upload()
        {
            var bytes = await Request.Content.ReadAsByteArrayAsync();
            if (bytes == null || bytes.Length == 0) return BadRequest("Empty body.");

            var tmp = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString("N") + ".dwxmz");
            File.WriteAllBytes(tmp, bytes);
            try
            {
                var id = FlowsheetStore.Create(auto => auto.LoadFlowsheet2(tmp));
                return Ok(new { id });
            }
            finally { TryDelete(tmp); }
        }

        // GET /api/flowsheet/{id}/objects  -> all simulation objects
        [HttpGet, Route("{id}/objects")]
        public IHttpActionResult ListObjects(string id)
        {
            var s = FlowsheetStore.Get(id);
            if (s == null) return NotFound();

            var objs = s.Flowsheet.SimulationObjects.Values.Select(o => new ObjectSummary
            {
                Name = o.Name,
                Tag = o.GraphicObject != null ? o.GraphicObject.Tag : o.Name,
                Type = o.GetType().Name,
                Calculated = o.Calculated,
                ErrorMessage = o.ErrorMessage
            }).ToList();

            return Ok(objs);
        }

        // GET /api/flowsheet/{id}/objects/{tag}  -> every property + current value
        [HttpGet, Route("{id}/objects/{tag}")]
        public IHttpActionResult GetObject(string id, string tag)
        {
            var s = FlowsheetStore.Get(id);
            if (s == null) return NotFound();

            var obj = s.Flowsheet.GetFlowsheetSimulationObject(tag);
            if (obj == null) return NotFound();

            var props = obj.GetProperties(PropertyType.ALL)
                .Select(p => new PropertyValue { Property = p, Value = SafeGet(obj, p) })
                .ToList();

            return Ok(new { tag, calculated = obj.Calculated, error = obj.ErrorMessage, properties = props });
        }

        // POST /api/flowsheet/{id}/property  -> set one property on one object
        [HttpPost, Route("{id}/property")]
        public IHttpActionResult SetProperty(string id, [FromBody] SetPropertyRequest req)
        {
            var s = FlowsheetStore.Get(id);
            if (s == null) return NotFound();
            if (req == null || string.IsNullOrEmpty(req.Tag) || string.IsNullOrEmpty(req.Property))
                return BadRequest("tag and property are required.");

            var obj = s.Flowsheet.GetFlowsheetSimulationObject(req.Tag);
            if (obj == null) return NotFound();

            lock (s.Lock)
            {
                var ok = obj.SetPropertyValue(req.Property, Coerce(req.Value));
                if (!ok) return BadRequest("SetPropertyValue rejected the value.");
            }
            return Ok();
        }

        // POST /api/flowsheet/{id}/solve  -> run the solver
        [HttpPost, Route("{id}/solve")]
        public IHttpActionResult Solve(string id)
        {
            var s = FlowsheetStore.Get(id);
            if (s == null) return NotFound();

            List<Exception> errors;
            lock (s.Lock) // serialize: engine is not thread-safe per flowsheet
            {
                errors = FlowsheetStore.Engine.CalculateFlowsheet4(s.Flowsheet);
            }

            return Ok(new SolveResult
            {
                Success = errors == null || errors.Count == 0,
                Errors = errors == null ? new List<string>() : errors.Select(e => e.Message).ToList()
            });
        }

        // GET /api/flowsheet/{id}/download  -> save and return the .dwxmz
        [HttpGet, Route("{id}/download")]
        public HttpResponseMessage Download(string id)
        {
            var s = FlowsheetStore.Get(id);
            if (s == null) return Request.CreateResponse(HttpStatusCode.NotFound);

            var tmp = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString("N") + ".dwxmz");
            try
            {
                lock (s.Lock) { FlowsheetStore.Engine.SaveFlowsheet2(s.Flowsheet, tmp); }
                var resp = Request.CreateResponse(HttpStatusCode.OK);
                resp.Content = new ByteArrayContent(File.ReadAllBytes(tmp));
                resp.Content.Headers.ContentType = new System.Net.Http.Headers.MediaTypeHeaderValue("application/octet-stream");
                resp.Content.Headers.ContentDisposition =
                    new System.Net.Http.Headers.ContentDispositionHeaderValue("attachment") { FileName = id + ".dwxmz" };
                return resp;
            }
            finally { TryDelete(tmp); }
        }

        [HttpDelete, Route("{id}")]
        public IHttpActionResult Delete(string id)
        {
            return FlowsheetStore.Remove(id) ? (IHttpActionResult)Ok() : NotFound();
        }

        // --- helpers ---

        private static object SafeGet(ISimulationObject obj, string prop)
        {
            try { return obj.GetPropertyValue(prop); }
            catch { return null; }
        }

        // JSON numbers arrive as long/double; DWSIM numeric props want Double.
        private static object Coerce(object v)
        {
            if (v == null) return null;
            if (v is long || v is int || v is float) return Convert.ToDouble(v);
            return v;
        }

        private static void TryDelete(string path)
        {
            try { if (File.Exists(path)) File.Delete(path); } catch { }
        }
    }
}
