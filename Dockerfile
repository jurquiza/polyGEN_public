FROM ubuntu:18.04
WORKDIR /scripts
RUN apt-get upgrade -y
RUN apt-get update -y
RUN apt-get install python3.6 -y
RUN apt-get install python3-pip -y #!
RUN apt-get update -y
RUN apt-get install git -y #!
RUN pip3 install jupyterlab==2.2.0
RUN pip3 install matplotlib==3.2.2
RUN pip3 install biopython==1.77
RUN pip3 install pysbol==2.3.3.post8
RUN pip3 install flask==1.1.2
RUN pip3 install flask-bootstrap==3.3.7.1
RUN pip3 install flask-wtf==0.14.3
RUN pip3 install git+https://github.com/Edinburgh-Genome-Foundry/icebreaker.git@e85c4f0121f4c8429eaf6ca3e2bd62a40dfb004e #using commit hash
WORKDIR /scripts
ENV FLASK_APP=/scripts/polygen.py
ENV FLASK_RUN_HOST=0.0.0.0
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
EXPOSE 5000
CMD ["flask","run"]
