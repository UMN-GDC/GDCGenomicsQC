# syntax=docker/dockerfile:1

#Python base image
FROM python:latest 

# Set working direcotry (also creates directory)
WORKDIR ~/app

# Install GCTA and plink
RUN curl -L https://zzz.bwh.harvard.edu/plink/dist/plink-1.08-x86_64.zip ; unzip plink-1.08-x86_64.zip ; rm plink-1.08-x86_64.zip
RUN curl -L https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip -o gcta.zip ; unzip gcta.zip ; rm gcta.zip

# Install ancestry prediction
git clone https://github.com/daviddaiweizhang/fraposa.git
cd fraposa
wget https://upenn.app.box.com/v/fraposa-demo/file/1170205195226

# Copy over test data
RUN mkdir Estimate
RUN mkdir Simulate

# copy files over to image
COPY Estimate Estimate 
COPY Simulate Simulate 

#load the python script and tell docker to run that script
#when someone tries to execute the container
ENTRYPOINT ["python3", "Estimate.py"]
