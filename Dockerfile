# Multi-stage Dockerfile for SAS Scatter2 Console Application with MCP Server
# Based on the docker commands from the documentation

# Stage 1: Build stage for the Qt application
FROM g76r/qt6-builder:qt-6.8.3-debug AS builder

# Set working directory
WORKDIR /home/user/project

# Copy source code into the container
COPY Source/ .

# Build the application
RUN rm -rf build && \
    mkdir build && \
    cd build && \
    qmake ../sas_scatter2Cons.pro && \
    make

# Stage 2: MCP Server build stage
FROM ghcr.io/astral-sh/uv:python3.13-bookworm AS mcp-builder

# Set working directory
WORKDIR /app/crystal_scatter_mcp

# Copy MCP server source code
COPY crystal_scatter_mcp/ .

# Install dependencies using uv
RUN uv sync --frozen

# Stage 3: Runtime stage
FROM g76r/qt6-runner:qt-6.8.3-debug AS runtime

# Update package database and install required libraries including uv
RUN apt-get update && \
    apt-get install -y \
    libxkbcommon-dev \
    libgl-dev \
    python3 \
    python3-venv \
    curl \
    && curl -LsSf https://astral.sh/uv/install.sh | sh \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the compiled executable from build stage
COPY --from=builder /home/user/project/build/sas_scatter2Cons /bin/sas_scatter2Cons

# Make the executable runnable
RUN chmod +x /bin/sas_scatter2Cons

# Copy MCP server files from mcp-builder stage
COPY --from=mcp-builder /app/crystal_scatter_mcp /app/crystal_scatter_mcp

# Set working directory for runtime
WORKDIR /app

# Create a startup script for the MCP server
COPY <<EOF /app/start-mcp.sh
#!/bin/bash
set -e

# ASCII art banner
cat << 'BANNER'
 ____              _   _             __  __  ____ ____  
/ ___|  __ _  ___ | |_| |_ ___ _ __ |  \/  |/ ___|  _ \ 
\___ \ / _` |/ _ \| __| __/ _ \ '__|| |\/| | |   | |_) |
 ___) | (_| | (_) | |_| ||  __/ |   | |  | | |___|  __/ 
|____/ \__,_|\___/ \__|\__\___|_|   |_|  |_|\____|_|    
                                                        
Crystal Scatter MCP Server                              
BANNER

echo "🚀 Starting MCP server..."
echo "📍 Working directory: /app/crystal_scatter_mcp"
echo "🔧 SAS Scatter2 Console: /bin/sas_scatter2Cons"
echo "🌐 Server will be available on port 8000"
echo ""

# Change to MCP server directory
cd /app/crystal_scatter_mcp

# Start the MCP server
exec uv run python server/mcp_server.py
EOF

RUN chmod +x /app/start-mcp.sh

# Expose port for MCP server (adjust if needed)
EXPOSE 8000

# Default command runs the MCP server
CMD ["/app/start-mcp.sh"]