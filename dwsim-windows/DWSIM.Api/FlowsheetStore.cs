using System;
using System.Collections.Concurrent;
using DWSIM.Automation;
using DWSIM.Interfaces;

namespace DWSIM.Api
{
    // One shared engine; one IFlowsheet per session, each with its own lock.
    // ponytail: in-memory, process-lifetime store. Swap for Redis/DB when you need >1 instance or persistence.
    public static class FlowsheetStore
    {
        private static Automation3 _auto;
        private static readonly ConcurrentDictionary<string, Session> _sessions =
            new ConcurrentDictionary<string, Session>();

        public class Session
        {
            public IFlowsheet Flowsheet;
            public readonly object Lock = new object(); // engine is NOT thread-safe per flowsheet
        }

        public static void Init()
        {
            _auto = new Automation3(); // loads property packages — expensive, once
        }

        public static string Create(Func<Automation3, IFlowsheet> build)
        {
            var id = Guid.NewGuid().ToString("N");
            _sessions[id] = new Session { Flowsheet = build(_auto) };
            return id;
        }

        public static Session Get(string id)
        {
            Session s;
            return _sessions.TryGetValue(id, out s) ? s : null;
        }

        public static bool Remove(string id)
        {
            Session s;
            return _sessions.TryRemove(id, out s);
        }

        public static Automation3 Engine => _auto;
    }
}
