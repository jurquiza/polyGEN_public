FROM ubuntu:20.04
WORKDIR /scripts
RUN apt-get upgrade -y
RUN apt-get update -y

RUN apt-get install -y --no-install-recommends software-properties-common \
    libsm6 libxext6 libxrender-dev curl \
    && rm -rf /var/lib/apt/lists/*

RUN echo "*** Installing Python v3.8***" && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get install -y build-essential python3.8 python 3.8-dev

RUN apt-get install python3-pip -y #!
RUN apt-get update -y
RUN apt-get install git -y #!

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

ENV FLASK_APP=/scripts/polygen.py
ENV FLASK_RUN_HOST=0.0.0.0
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
EXPOSE 5000
CMD ["flask","run"]
