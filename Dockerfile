# Multi-stage Dockerfile: Crystal Scatter CLI + MCP Server
FROM g76r/qt6-builder:qt-6.8.3-debug AS sas_scatter-builder

# Build SAS Scatter CLI from Source directory
WORKDIR /sas_scatter
COPY Source/ /sas_scatter/

# Clean build
RUN rm -rf build && \
    mkdir build && \
    cd build && \
    qmake ../sas_scatter2Cons.pro && \
    make -j$(nproc) && \
    strip sas_scatter2Cons

# Final runtime stage - use official uv image
FROM ghcr.io/astral-sh/uv:python3.13-alpine

# Install system dependencies (Qt runtime + utilities)
RUN apk add --no-cache \
    qt6-qtbase \
    curl \
    ca-certificates

# Create app directory and user
RUN addgroup -S mcpuser && adduser -S mcpuser -G mcpuser -u 1000
WORKDIR /app

# Copy SAS Scatter executable from builder
COPY --from=sas_scatter-builder /sas_scatter/build/sas_scatter2Cons /usr/local/bin/sas_scatter2Cons
RUN chmod +x /usr/local/bin/sas_scatter2Cons

# Copy MCP server source code
COPY crystal_scatter_mcp/ ./

# Install the package and dependencies
RUN uv pip install --system -e .

# Create necessary directories
RUN mkdir -p data cache logs && \
    chown -R mcpuser:mcpuser /app

# Switch to non-root user
USER mcpuser

# Expose MCP server port
EXPOSE 8000

# Health check endpoint
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/mcp || exit 1

# Default command
CMD ["python", "-m", "server.mcp_server", "http", "--host", "0.0.0.0", "--port", "8000"]