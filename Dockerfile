# Use a standard Python image
FROM python:3.11-slim

# Install system dependencies required by spatial libraries (e.g. OpenCV, gcc)
RUN apt-get update && apt-get install -y \
    build-essential \
    libgl1 \
    libglib2.0-0 \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# 1. Install EXACT versions of pip and setuptools FIRST
RUN pip install --upgrade pip==26.0.1 setuptools==81.0.0 wheel

# 2. Install MuSpan 
# Define variables that we will pass during the build command
ARG MUSPAN_USER
ARG MUSPAN_PASS

# Run the installation by injecting the credentials into the URL
RUN curl -L -u "${MUSPAN_USER}:${MUSPAN_PASS}" -o muspan.zip "https://docs.muspan.co.uk/code/latest.zip" \
    && pip install muspan.zip \
    && rm muspan.zip

# Copy requirement list 
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy all python files
COPY . /app

# The default command will run the pipeline. 
ENTRYPOINT ["python", "modules/__main__.py"]
CMD ["config_CosMx.toml"]