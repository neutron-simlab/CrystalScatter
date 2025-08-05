"""Main MCP server application with HTTP transport"""
import typer
import subprocess
from fastmcp import FastMCP
from typing import Dict, Any
from pathlib import Path
from pydantic import BaseModel
from typing import Dict, Any, Optional
from pathlib import Path

class CrystalScatterCLI:
    """Handles Crystal Scatter CLI application execution with proper error handling"""
    
    def __init__(self, executable_name: str = "sas_scatter2Cons"):
        self.executable_name = executable_name
        
    def check_executable_exists(self) -> bool:
        """Check if the CLI executable is available"""
        try:
            result = subprocess.run(
                ['which', self.executable_name], 
                capture_output=True, 
                text=True
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def get_executable_path(self) -> Optional[str]:
        """Get the full path to the executable"""
        try:
            result = subprocess.run(
                ['which', self.executable_name], 
                capture_output=True, 
                text=True
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return None
        except Exception:
            return None
    
    def execute_with_str(
        self, 
        kvstr: str, 
        timeout: int = 30,
    ) -> Dict[str, Any]:
        """
        Execute Crystal Scatter CLI app with key value pairs string and return structured output
        
        Args:
            kvstr: Key value pair string delimited by ; (semicolon)
            timeout: Maximum execution time in seconds
            
        Returns:
            Dictionary with execution results
        """
        
        try:
            # Build command
            cmd = [self.executable_name, "--mcpval", kvstr, "--mcpimg64", "--threads", "0" ]
            
            # Execute command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False
            )
            
             
            
            return {
                "success": result.returncode == 0,
                "exit_code": result.returncode,
                "base64_img": result.stdout.strip(),
                "stderr": result.stderr,
                "command": ' '.join(cmd),
            }
            
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "exit_code": -1,
                "stdout": "",
                "stderr": f"Command timed out after {timeout} seconds",
                "command": ' '.join(cmd) if 'cmd' in locals() else "unknown", # type: ignore
            }
        except Exception as e:
            return {
                "success": False,
                "exit_code": -1,
                "stdout": "",
                "stderr": f"Execution error: {str(e)}",
                "command": ' '.join(cmd) if 'cmd' in locals() else "unknown", # type: ignore
            }
        
    def execute_with_file(
        self, 
        input_path: Path,
        output_path: Path = Path.home(),
        timeout: int = 30,
    ) -> Dict[str, Any]: 
        """
        Execute Crystal Scatter CLI app with an input file and return structured output as well output image file.
    
        Args:
            input_path: A path to the crystal scatter input file
            output_path: A path where the image from the simulaiton should be stored
            timeout: Maximum execution time in seconds
            
        Returns:
            Dictionary with execution results
        """
        
        try:
            # Build command
            cmd = [self.executable_name, "--mcpinp", str(input_path), "--mcpimg",  str(output_path), "--threads", "0" ]
            
            # Execute command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False
            )
            
            return {
                "success": result.returncode == 0,
                "exit_code": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": ' '.join(cmd),
            }
            
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "exit_code": -1,
                "stdout": "",
                "stderr": f"Command timed out after {timeout} seconds",
                "command": ' '.join(cmd) if 'cmd' in locals() else "unknown", # type: ignore
            }
        except Exception as e:
            return {
                "success": False,
                "exit_code": -1,
                "stdout": "",
                "stderr": f"Execution error: {str(e)}",
                "command": ' '.join(cmd) if 'cmd' in locals() else "unknown", # type: ignore
            }

class CLIBase64Image(BaseModel):
    success:bool
    exit_code:int
    base64_image:str
    stderr: str
    command: str

# Initialize FastMCP server with HTTP transport
crystal_scatter_mcp = FastMCP("Remote Crystal Structure CLI MCP Server")
crystal_scatter_cli = CrystalScatterCLI()

# MCP Tools
@crystal_scatter_mcp.tool()
def run_cli_app_with_str(kvstr: str, timeout: int = 30) -> Dict[str, Any]:
    """
    Execute the Crystal Scatter CLI application with key value pairs string and capture output.
    
    Args:
        kvstr: String containing key value pairs of crystal scatter input parameters
        timeout: Maximum execution time in seconds
    
    Returns:
        String containing the captured output
    """
    result = crystal_scatter_cli.execute_with_str(kvstr, timeout)

    valid_result = CLIBase64Image(
        success=result['success'],
        exit_code=result['exit_code'],
        base64_image=result['base64_img'],
        stderr=result['stderr'],
        command=result['command'],
        )
    
    return valid_result.model_dump()

@crystal_scatter_mcp.tool()
def run_cli_app_with_file(input_path: Path, output_path: Path = Path.home(), timeout: int = 30) -> str:
    """
    Execute Crystal Scatter CLI app with an input file and return structured output as well output image file.
        
    Args:
        input_path: A path to the crystal scatter input file
        output_path: A path where the image from the simulaiton should be stored    
        timeout: Maximum execution time in seconds
    
    Returns:
        String containing the captured output
    """
    result = crystal_scatter_cli.execute_with_file(input_path, output_path, timeout)
    
    # Format output for MCP
    output_parts = []
    if result["stdout"]:
        output_parts.append(f"Output:\n{result['stdout']}")
    if result["stderr"]:
        output_parts.append(f"Errors:\n{result['stderr']}")
    output_parts.append(f"Exit code: {result['exit_code']}")
    output_parts.append(f"Command: {result['command']}")
    
    return "\n\n".join(output_parts)

# Internal health check function (not an MCP tool)
def _internal_health_check() -> str:
    """Internal health check function"""
    executable_found = crystal_scatter_cli.check_executable_exists()
    executable_path = crystal_scatter_cli.get_executable_path()
    
    if executable_found and executable_path:
        return f"✅ Server healthy. CLI app found at: {executable_path}"
    else:
        return "❌ CLI app not found in PATH"

@crystal_scatter_mcp.tool()
def health_check() -> str:
    """Check if the server and CLI app are working"""
    return _internal_health_check()

@crystal_scatter_mcp.tool()
def get_server_info() -> str:
    """Get information about the MCP server"""
    return "Remote Crystal Structure CLI MCP Server - HTTP Transport Mode"

# CLI interface
def main():
    """Main entry point with CLI interface"""
    app_cli = typer.Typer()
    
    @app_cli.command()
    def http(
        port: int = typer.Option(8000, help="Port to run HTTP server on"),
        host: str = typer.Option("0.0.0.0", help="Host to bind to")
    ):
        """Run server in HTTP transport mode"""
        print(f"🚀 Starting MCP server with HTTP transport on {host}:{port}")
        
        # Configure FastMCP to use HTTP transport
        crystal_scatter_mcp.run(transport="http", host=host, port=port)
    

    @app_cli.command()
    def health():
        """Check server health"""
        result = _internal_health_check()
        print(result)
    
    app_cli()

if __name__ == "__main__":
    main()