FROM kbase/sdkbase2:python
MAINTAINER Alex Brace
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

# update security libraries in the base image
RUN pip install --upgrade pip setuptools wheel cffi && \
    pip install --upgrade pyopenssl ndg-httpsclient && \
    pip install --upgrade pyasn1 requests 'requests[security]' && \
    pip install coverage networkx cython

# Install forked version of optlang and cobrapy to add
# additional solver support COINOR-CBC,CLP and OSQP
RUN mkdir deps && cd deps && \
    git clone https://github.com/braceal/optlang.git && \
    cd optlang && git checkout test/coinor-cbc_osqp && cd .. && \
    git clone https://github.com/braceal/cobrapy.git && \
    cd cobrapy && git checkout feature/coinor-cbc_osqp && cd .. && \
    pip install optlang/ && \
    pip install cobrapy/ --ignore-installed && cd .. && \
    pip install cobrakbase==0.2.7 --ignore-installed && \
    pip install Jinja2 

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
