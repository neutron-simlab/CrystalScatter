# Crystal Scatter 

## Overview

This project combines the **Crystal Scatter CLI application** with a **Model Context Protocol (MCP) server**, enabling LLMs and other MCP clients to execute crystal structure calculations remotely through a simple HTTP interface.

```
LLM-backed Chatbot + MCP Client → HTTP (MCP Protocol) → FastMCP Server → Crystal Scatter CLI
```

## Crystal Scatter Background

Crystal Scatter is a high-performance small-angle scattering simulation software developed at Forschungszentrum Jülich.

- **Original Pascal version**: Stephan Förster (s.foerster@fz-juelich.de)
- **C++ with CUDA conversion**: Michael Wagener (m.wagener@fz-juelich.de)
- **Model Context Protocol (MCP) Server Integration**: Ahmad Zainul Ihsan (a.ihsan@fz-juelich.de)
- **License**: Open Source © (2023-2025) Forschungszentrum Jülich GmbH, JCNS-1
- **Scientific Publication**: [Nature Scientific Reports](https://www.nature.com/articles/s41598-023-27558-8)

### Features
- High-performance crystal structure simulations
- CUDA GPU acceleration (Linux only)
- Multi-threaded CPU processing
- Support for various particle geometries (spheres, disks, cylinders)
- Multiple scattering models (homogeneous, core-shell)

### Used Libraries
The given version numbers are used during development and are working. Other library versions may work too. To configure the used path of the external libraries edit the `sas_scatter2.pri` file in the source directory.

#### Mandatory
- Qt 5.14.x (first development)
- Qt 6.7.2 (current version)

#### Optional / Windows
- Qwt-6.3.0
- fftw-3.3.5-dll64
- hdf5-1.12.1 (CMake build)

#### Optional / Linux
- cuda-rhel8-11-0-local
- qwt-qt5-devel.x86_64 (6.1.5-5.el8)
- fftw-devel.x86_64 (3.3.5-11.el8)
- fftw-libs.x86_64 (3.3.5-11.el8)
- hdf-devel.x86_64 (4.2.14-5.el8)
- or download and install manually from hdf5 website (checked with 1.12.1 and 1.14.4-2)

**Note**: If you want to use a CUDA version newer than 11.7, you have to set the environment variable `CUDA_MODULE_LOADING=EAGER` before compile and launch to disable the lazy loading feature.

## MCP Server Features

- **HTTP Transport**: Uses MCP protocol over HTTP for easy deployment and scaling
- **Dockerized Deployment**: Single container with both Crystal Scatter CLI and MCP server
- **Health Monitoring**: Built-in health checks and status monitoring
- **Async Support**: Full async/await support for concurrent operations
- **LLM Integration**: Ready for integration with AI chatbots and assistants

## Quick Start with Docker

### Prerequisites
- Docker and Docker Compose
- Python 3.13+ (for running the client)
- uv or pip (for client dependencies)

### Setup

1. **Repository structure:**
```bash
crystal-scatter/
├── Source/                     # Crystal Scatter source files
│   ├── sas_scatter2Cons.pro   # Qt project file
│   └── (other source files)
├── crystal_scatter_mcp/       # MCP server integration
│   ├── client/                # MCP client code
│   ├── server/                # MCP server code
│   ├── pyproject.toml
│   └── uv.lock
├── Dockerfile                  # Container definition
├── docker-compose.yml         # Service orchestration
```

2. **Build and start the MCP server:**
```bash
# From the repository root
docker-compose up -d

# Check if server is running
docker-compose ps

# Check server health
curl http://localhost:8000/mcp
```

3. **Install client dependencies (for testing):**

**Using uv (recommended):**
```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create and activate virtual environment
uv venv
source .venv/bin/activate  # Linux/Mac
# or .venv\Scripts\activate  # Windows

# Install client dependencies
uv add fastapi fastmcp pydantic typer uvicorn
```

**Using pip:**
```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# or venv\Scripts\activate  # Windows

# Install client dependencies
pip install fastapi>=0.116.0 fastmcp>=2.10.2 pydantic>=2.11.7 typer>=0.16.0 uvicorn>=0.35.0
```

4. **Test the MCP server:**
```bash
# Run quick test (from repository root)
python crystal_scatter_mcp/client/mcp_client.py quick-test

# Or just health check
python crystal_scatter_mcp/client/mcp_client.py health
```

## Architecture

### Container Approach
The Docker container uses a multi-stage build:

1. **Build Stage**: Compiles Crystal Scatter CLI using Qt6 builder with additional dependencies (FFTW, HDF5)
2. **Runtime Stage**: Python environment with manually installed dependencies + compiled CLI

### Available MCP Tools

1. **`run_cli_app_with_str`**
   - Execute Crystal Scatter CLI with parameter string
   - Parameters: `kvstr` (string), `timeout` (int, default: 30)
   - Returns: Base64 encoded simulation image and metadata

2. **`run_cli_app_with_file`**
   - Execute Crystal Scatter CLI with input file
   - Parameters: `input_path`, `output_path`, `timeout`
   - Returns: Simulation results and output file path

3. **`health_check`**
   - Check server and CLI application status
   - Returns: Health status message

4. **`get_server_info`**
   - Get server information
   - Returns: Server description

## Usage Examples

### Direct MCP Client Usage

**Prerequisites**: Install client dependencies (see setup section above)

```bash
# From repository root
python crystal_scatter_mcp/client/mcp_client.py quick-test
python crystal_scatter_mcp/client/mcp_client.py health
```

**Or use the provided test client:**
```bash
# From repository root
python crystal_scatter_mcp/client/mcp_client.py quick-test
python crystal_scatter_mcp/client/mcp_client.py health
```

### LLM Integration
The MCP server is designed to work with LLM chatbots that support tool calling:

```python
# Example with OpenAI
tools = [
    {
        "type": "function",
        "function": {
            "name": "run_crystal_scatter_simulation",
            "description": "Run crystal scatter simulation",
            "parameters": {
                "type": "object",
                "properties": {
                    "parameters": {"type": "string", "description": "Simulation parameters"}
                }
            }
        }
    }
]

# LLM can now call crystal scatter simulations based on natural language
```

```python
# Example with Langchain + Langchain MCP adapter
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_openai import ChatOpenAI

client = MultiServerMCPClient(
    {
        "crystal_scatter_cli_tools": {
            "command": "python",
            # Make sure to update to the full absolute path to your math_server.py file
            "args": ["/crystal_scatter_mcp/server/mcp_server.py", "http"],
            "transport": "http",
)
tools = await client.get_tools()
llm = ChatOpenAI(model=...).bind_tools(tools)
```



## Configuration

### Server Configuration
The server runs inside Docker container with these settings:
- **Port**: 8000 (exposed to host)
- **Dependencies**: Installed via pip (fastapi, fastmcp, pydantic, typer, uvicorn)
- **Executable**: Crystal Scatter CLI available at `/bin/sas_scatter2Cons`

### Client Configuration
The client runs on the host machine and requires:
- **Python**: 3.13+ (minimum required version)
- **Dependencies**: fastapi, fastmcp, pydantic, typer, uvicorn
- **Server URL**: `http://localhost:8000/mcp`

### Environment Variables
```bash
# Optional: Create .env file for custom settings
MCP_SERVER_HOST=0.0.0.0
MCP_SERVER_PORT=8000
SIMULATION_TIMEOUT=60
```

## Crystal Scatter Dependencies

The container includes all necessary dependencies:

### Runtime Dependencies (Included in Docker)
- Qt 6.8.3 (runtime libraries)
- Python 3.13+
- FastMCP server
- FFTW (Fast Fourier Transform library)
- HDF5 (Hierarchical Data Format)
- System utilities (curl, etc.)

### Build Dependencies (Build-time only)
- Qt 6.8.3 development environment
- GCC/G++ compiler
- CMake and build tools
- pkg-config

## Integration Examples

### Chatbot Integration
The MCP server is designed to work with various LLM chatbots:

- **OpenAI GPT**: Using function calling
- **Anthropic Claude**: Using tool use
- **Langchain**: ChatModel bind tools
- **Custom chatbots**: Any MCP-compatible client

### API Endpoints
- **MCP Protocol**: `http://localhost:8000/mcp`
  
## Documentation

- **Crystal Scatter Documentation**: [Documentation](./Documentation) (up to date)
- **Scientific Background**: [HyperSeries_SI_rev.pdf](./HyperSeries_SI_rev.pdf)
- **MCP Specification**: [Model Context Protocol](https://modelcontextprotocol.io/)
- **FastMCP Documentation**: [FastMCP GitHub](https://github.com/jlowin/fastmcp)

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test with the included demo client
5. Submit a pull request

## License

- **Crystal Scatter**: Open Source © (2023-2025) Forschungszentrum Jülich GmbH, JCNS-1

## Links

- [Crystal Scatter CLI Documentation](https://github.com/neutron-simlab/CrystalScatter/Documentation)
- [FastMCP Documentation](https://github.com/jlowin/fastmcp)
- [Model Context Protocol Specification](https://modelcontextprotocol.io/)
- [Forschungszentrum Jülich](https://www.fz-juelich.de/)
- [Scientific Publication](https://www.nature.com/articles/s41598-023-27558-8)

## Quick Commands Reference

```bash
# Build and start server
docker-compose up -d

# Install client dependencies (choose one)
uv venv && source .venv/bin/activate && uv add fastapi fastmcp pydantic typer uvicorn
# OR
python3 -m venv venv && source venv/bin/activate && pip install fastapi fastmcp pydantic typer uvicorn

# Test the MCP server
python crystal_scatter_mcp/client/mcp_client.py quick-test

# Check health
python crystal_scatter_mcp/client/mcp_client.py health

# View logs
docker-compose logs -f crystal-scatter-mcp

# Stop services
docker-compose down
```

## Common Workflow

1. **Start the server**: `docker-compose up -d`
2. **Set up client environment**: Install Python dependencies 
3. **Run simulations**: Use the test client or write your own MCP client
4. **Monitor**: Check logs and health endpoints
5. **Integrate**: Connect your LLM chatbot to `http://localhost:8000/mcp`
