# Start from a clean, stable Linux base with Miniconda installed
FROM continuumio/miniconda3:latest

# Set a working directory inside the container
WORKDIR /app

# Copy the exported environment file into the container image
COPY environment.yml /app/

# Install the Conda environment from the YAML file
# -n my_app_env names the environment (optional, but good practice)
RUN conda install -y mamba -c conda-forge && \
    mamba env create -f environment.yml -n my_app_env && \
    # Clean up unnecessary files to reduce the final image size
    conda clean --all -f -y

# Set the PATH to include the new environment's bin directory
ENV PATH="/opt/conda/envs/my_app_env/bin:$PATH"

# Copy your local source code, scripts, or configuration files
# Assuming your main script is in the local 'src' directory
COPY src/ /app/src/
COPY data/ /app/src/data/
COPY results/ /app/src/results/
COPY Run.sh /app/

ENTRYPOINT= ["/app/Run.sh"]


# Define the command that runs when the container is executed
# This activates the environment and runs your main script
CMD ["--help"]

