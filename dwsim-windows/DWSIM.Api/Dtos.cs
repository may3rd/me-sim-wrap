using System.Collections.Generic;

namespace DWSIM.Api
{
    public class SetPropertyRequest
    {
        public string Tag { get; set; }      // object tag, e.g. "Feed"
        public string Property { get; set; } // DWSIM property id, e.g. "PROP_MS_0" (temperature)
        public object Value { get; set; }
    }

    public class ObjectSummary
    {
        public string Name { get; set; }     // internal id (dictionary key)
        public string Tag { get; set; }      // user-facing tag
        public string Type { get; set; }
        public bool Calculated { get; set; }
        public string ErrorMessage { get; set; }
    }

    public class PropertyValue
    {
        public string Property { get; set; }
        public object Value { get; set; }
    }

    public class SolveResult
    {
        public bool Success { get; set; }
        public List<string> Errors { get; set; } = new List<string>();
    }
}
