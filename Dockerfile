# Stage 1: Build Qt app
FROM g76r/qt6-builder:qt-6.8.3-debug AS builder

USER root
RUN rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y \
    libfftw3-dev \
    libhdf5-dev \
    pkg-config \
    ca-certificates \
    gnupg

WORKDIR /home/user/project
COPY Source/ .

RUN mkdir build && \
    cd build && \
    qmake ../sas_scatter2Cons.pro && \
    make

# Stage 2: Runtime
FROM g76r/qt6-runner:qt-6.8.3-debug AS runtime

USER root
RUN apt-get update && \
    apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    libxkbcommon-dev \
    libgl-dev \
    curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy the Qt console app
COPY --from=builder /home/user/project/build/sas_scatter2Cons /bin/sas_scatter2Cons
RUN chmod +x /bin/sas_scatter2Cons

# App setup
WORKDIR /app/crystal_scatter_mcp
COPY crystal_scatter_mcp/ .

# Create a venv and install dependencies manually from pyproject
RUN python3 -m venv /app/venv && \
    /app/venv/bin/pip install --upgrade pip && \
    /app/venv/bin/pip install fastapi>=0.116.0 fastmcp>=2.10.2 pydantic>=2.11.7 typer>=0.16.0 uvicorn>=0.35.0

# Expose port
EXPOSE 8000

# Run directly
CMD ["/app/venv/bin/python", "server/mcp_server.py", "http"]
