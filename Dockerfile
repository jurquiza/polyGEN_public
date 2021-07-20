FROM ubuntu:latest
WORKDIR /scripts
RUN apt-get upgrade -y
RUN apt-get update -y
RUN apt-get install python3 -y
RUN apt-get install python3-pip -y #!
RUN apt-get update -y
RUN apt-get install git -y #!
RUN pip3 install numpy==1.21.1
RUN pip3 install pandas==1.3.0
RUN pip3 install biopython==1.77
RUN pip3 install flask==1.1.2
RUN pip3 install flask-bootstrap==3.3.7.1
RUN pip3 install flask-wtf==0.14.3
WORKDIR /scripts
ENV FLASK_APP=/scripts/polygen.py
ENV FLASK_RUN_HOST=0.0.0.0
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
EXPOSE 5000
CMD ["flask","run"]
