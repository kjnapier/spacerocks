# Dockerfile courtesy of JJ Kavelaars

FROM ubuntu:latest

USER root
RUN apt-get update  -y -q
# canfar science portal bits.
RUN apt-get install -q -y sssd libnss-sss libpam-sss
# system settings and permissions
ADD skaha/nsswitch.conf /etc/
ADD skaha/nofiles.conf /etc/security/limits.d/
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf
# generate missing dbus uuid (issue #47)
## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN dbus-uuidgen --ensure

# application 
RUN apt install -q -y build-essential libssl-dev libffi-dev python3-dev
RUN apt install -q -y python3-pip
RUN apt install -q -y python3-numpy
RUN apt install -q -y python3-scipy
RUN pip3 install vos cadcdata cadctap cadcutils

# RUN apt-get install -y g++
# RUN apt-get install -y libpcre3 libpcre3-dev
RUN apt install -y swig

ARG BUILDDIR=/opt/spacerocks
RUN mkdir -p ${BUILDDIR}
WORKDIR ${BUILDDIR}
# COPY swig-4.0.2.tar.gz ./
# RUN tar -xzf swig-4.0.2.tar.gz
# WORKDIR ${swig-4.0.2}
# RUN ./configure && make && make install
# WORKDIR ${BUILDDIR}
COPY  ./ ./
RUN python3 -m pip install -r requirements.txt
RUN python3 setup.py install 
RUN pip install jupyter

CMD [ "jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root" ]

