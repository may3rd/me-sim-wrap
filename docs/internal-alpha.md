# Internal alpha deployment

The service is **private-network only**. It has no application authentication and must not be exposed through a public IP, public load balancer, or public ingress.

Build the image for the target architecture, attach it only to the internal application network, and allow requests only from the in-house web app. The image intentionally excludes `dwsim-windows/`; DWSIM remains a Windows reference oracle, not a runtime dependency.

Do not publish port 8000 on a public host. Add bearer-token authentication before any network exposure beyond the private application network.
