# Start from a clean, stable Linux base with Miniconda installed
FROM continuumio/miniconda3:latest

# Set a working directory inside the container
WORKDIR /app

# Copy the exported environment file into the container image
COPY environment.yml /app/

# Install the Conda environment from the YAML file
# -n my_app_env names the environment (optional, but good practice)
RUN conda install -y mamba -c conda-forge && \
    mamba env create -f environment.yml && \
    # Clean up unnecessary files to reduce the final image size
    conda clean --all -f -y


# Set the PATH to include the new environment's bin directory
ENV PATH="/opt/conda/envs/gdcPipeline/bin:$PATH"

# Install softwares not managed by conda
RUN pip3 install git+https://github.com/liguowang/CrossMap.git
RUN mkdir -p /opt/genotypeHarmonizer
RUN mkdir -p /opt/primus
RUN curl -SL -o https://github.com/molgenis/systemsgenetics/releases/download/GH_1.4.28/GenotypeHarmonizer-1.4.28-SNAPSHOT-dist.tar.gz | tar -xzC /opt/GenotypeHarmonizer
RUN curl -SL -o https://primus.gs.washington.edu/docroot/versions/PRIMUS_v1.9.0.tgzhttps://primus.gs.washington.edu/docroot/versions/PRIMUS_v1.9.0.tgz | tar -xzC /opt/primus --strip-components=1

# Wrapper script in /usr/bin
RUN echo '#!/bin/sh' > /usr/bin/GenotypeHarmonizer \
    && echo 'exec java -jar /opt/GenotypeHarmonizer/GenotypeHarmonizer.jar "$@"' >> /usr/bin/my-command \
    && chmod +x /usr/bin/
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

