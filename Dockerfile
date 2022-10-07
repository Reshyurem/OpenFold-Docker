FROM nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu18.04

# metainformation
LABEL org.opencontainers.image.version = "1.0.0"
LABEL org.opencontainers.image.authors = "Gustaf Ahdritz"
LABEL org.opencontainers.image.source = "https://github.com/aqlaboratory/openfold"
LABEL org.opencontainers.image.licenses = "Apache License 2.0"
LABEL org.opencontainers.image.base.name="docker.io/nvidia/cuda:10.2-cudnn8-runtime-ubuntu18.04"

RUN apt-key del 7fa2af80
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub

RUN apt-get update && apt-get install -y wget curl unzip libxml2 cuda-minimal-build-11-3 libcusparse-dev-11-3 libcublas-dev-11-3 libcusolver-dev-11-3 git hmmer vim

RUN wget -P /tmp \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /opt/openfold/environment.yml

# installing into the base environment since the docker container wont do anything other than run openfold
RUN conda env update -n base --file /opt/openfold/environment.yml && conda clean --all

RUN mkdir -m 777 --parents /tmp/ramdisk
# RUN mount --bind -t tmpfs -o size=9G ramdisk /tmp/ramdisk

RUN wget -q -P /content https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -qq awscliv2.zip
RUN ./aws/install
RUN rm -rf awscliv2.zip aws

RUN apt-get update && apt-get install -y python3 python3-pip
# RUN pip3 install --upgrade setuptools
# RUN pip3 install --upgrade pip

RUN rm -rf openfold
RUN git clone https://github.com/aqlaboratory/openfold openfold 2>&1 1> /dev/null
RUN mkdir -p /content/openfold/openfold/resources
RUN cp -f /content/stereo_chemical_props.txt /content/openfold/openfold/resources
RUN pip3 install -q ./openfold

# RUN ls -l /opt/conda/lib/python3.7/site-packages/

# RUN cd /opt/conda/lib/python3.7/site-packages/
# RUN patch -p0 < /openfold/lib/openmm.patch
# RUN cd -

COPY input.txt /openfold/input.txt
COPY openfold.py /openfold/main.py

RUN mkdir --parents /content/openfold/openfold/resources/openfold_params
RUN aws s3 cp --no-sign-request --region us-east-1 s3://openfold/openfold_params /content/openfold/openfold/resources/openfold_params --recursive

WORKDIR /openfold

# CMD python3 /openfold/openfold.py

# CMD uvicorn main:app --host 0.0.0.0 --port 8000