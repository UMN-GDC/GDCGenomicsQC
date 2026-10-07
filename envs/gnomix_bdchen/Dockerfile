FROM python:3.7-slim
RUN apt-get update && apt-get install -y gcc wget unzip bcftools
COPY requirements.txt .
ENV SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
RUN pip install -r requirements.txt
RUN mkdir /gnomix
RUN wget https://github.com/AI-sandbox/gnomix/archive/refs/heads/main.zip && \
    unzip main.zip && \
    mv gnomix-main/* /gnomix/ && \
    rm -rf main.zip /gnomix/demo && \
    cp /gnomix/config.yaml . && \
    cd /gnomix/