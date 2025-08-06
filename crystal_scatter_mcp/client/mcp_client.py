import asyncio
import typer
from fastmcp import Client
from pathlib import Path

async def run_test(host: str = "localhost", port: int = 8000):
    """Run quick health check and sample test on local MCP server"""
    
    print("🏠 Running local MCP quick test...")
    print(f"🌐 Connecting to: {host}:{port}/mcp")
    
    server_url = f"http://{host}:{port}/mcp"
    
    try:
        async with Client(server_url) as client:
            # Test basic connection
            print("🔄 Testing connection...")
            await client.ping()
            print("✅ Connection successful")
            
            # List available tools
            print("📋 Listing available tools...")
            tools = await client.list_tools()
            tool_names = [tool.name for tool in tools]
            print(f"🔧 Available tools: {', '.join(tool_names)}")
            
            # Health check
            print("🏥 Checking server health...")
            health_response = await client.call_tool("health_check", {})
            print(f"🏥 Health: {health_response.content[0].text}") # type: ignore
            
            # Prepare the key-value data
            keyval = [
                "Experiment_name=Disks",
                "base=0.0",
                "beamposx=57.0",
                "beamposy=31.0",
                "cbinterior=homogeneous",
                "cbparticle=disk",
                "ceff=0.01",
                "ceffcyl=0.0",
                "dbeta=0.4",
                "debyescherrer=False",
                "det=10.0",
                "gridpoints=100.0",
                "i0=1000.0",
                "iso=0.0",
                "length=2.0",
                "ltype=None",
                "ordis=isotropic",
                "peakpar=0.0",
                "phi=0.0",
                "pixelnox=128.0",
                "pixelnoy=128.0",
                "pixelx=1.0",
                "pixely=1.0",
                "qmax=2.0",
                "radius=4.0",
                "rbpara=False",
                "reff=0.0",
                "rotalpha=0.0",
                "rotphi=0.0",
                "rottheta=0.0",
                "sigma=0.1",
                "sigmal=0.1",
                "theta=0.0",
                "ucalpha=90.0",
                "ucbeta=90.0",
                "ucgamma=90.0",
                "ucn1=1.0",
                "ucn2=0.0",
                "ucn3=0.0",
                "ucpsi=0.0",
                "wavelength=0.154"
            ]

            input_file = Path("/home/ahmad/crystal-scatter-mcp/src/crystal_scatter_mcp/client/data/input.txt")
            output_file = Path("/home/ahmad/crystal-scatter-mcp/src/crystal_scatter_mcp/client/data/output.png")
            # Convert to single string
            kvstr = ';'.join(keyval)
            
            print(f"📤 Sending data: {len(keyval)} parameters")
            print(f"📝 Data preview: {kvstr[:100]}...")

            # Execute CLI command
            print("⚙️ Executing CLI command for the string input")
            response = await client.call_tool(
                "run_cli_app_with_str", 
                {
                    'kvstr': kvstr, 
                    'timeout': 30
                }
            )
            
            print("✅ Local test passed for string input")
            print("📋 Output base64 image:")
            print(response.data['base64_image']) # type: ignore


            # Execute CLI command using input and output file
            print("⚙️ Executing CLI command for the file input")
            response = await client.call_tool(
                "run_cli_app_with_file", 
                {
                    'input_path': input_file, 
                    'output_path': output_file
                }
            )
            
            print("✅ Local test passed for the file input")
            print("📋 Output:")
            print(response.content[0].text) # type: ignore
            
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        print(f"💡 Make sure the server is running: uv run python src/crystal_scatter_mcp/server/mcp_server.py http")

# CLI interface
app = typer.Typer()

@app.command()
def quick_test(
    host: str = typer.Option("0.0.0.0", help="Server hostname"),
    port: int = typer.Option(8000, help="Server port")
):
    """Run quick health check and sample test on local MCP server"""
    asyncio.run(run_test(host, port))

@app.command()
def health(
    host: str = typer.Option("0.0.0.0", help="Server hostname"),
    port: int = typer.Option(8000, help="Server port")
):
    """Just check server health"""
    async def health_check():
        server_url = f"http://{host}:{port}/mcp"
        try:
            async with Client(server_url) as client:
                await client.ping()
                health_response = await client.call_tool("health_check", {})
                print(f"🏥 Health: {health_response.content[0].text}") # type: ignore
        except Exception as e:
            print(f"❌ Health check failed: {str(e)}")
    
    asyncio.run(health_check())

def main():
    """Main entry point"""
    app()

if __name__ == "__main__":
    main()