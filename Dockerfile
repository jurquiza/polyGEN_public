FROM ubuntu:latest
RUN apt-get upgrade -y
RUN apt-get update -y
RUN apt-get install python3 -y
RUN apt-get install python3-pip -y
RUN apt-get update -y
RUN apt-get install git -y
RUN pip3 install jupyterlab
RUN pip3 install matplotlib
RUN pip3 install biopython
RUN pip3 install pysbol
RUN pip3 install flask
RUN pip3 install flask-bootstrap
RUN pip3 install flask-wtf
RUN pip3 install git+https://github.com/Edinburgh-Genome-Foundry/icebreaker.git
ADD ./commands/jn /bin/
RUN chmod +x /bin/jn
ADD ./commands/polygen /bin/
RUN chmod +x /bin/polygen
