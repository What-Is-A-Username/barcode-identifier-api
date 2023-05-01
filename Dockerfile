FROM python:3.8.10 AS pyvenv_image
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1

#TODO: postgresql-client apk needed?
#TODO: Remove .tmp-deps after install
RUN apt-get update && \
    apt-get install -y gcc \
    libc6-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    build-essential \
    python3-dev \
    linux-headers-amd64 \
    libpcre3-dev \
    libpq-dev \
    python3-dev \
    libffi-dev

# note: gcompat necessary to run our ncbi binaries on an alpine image
# RUN apk add --no-cache postgresql-dev \ 
#     gcc \
#     musl-dev \
#     curl-dev \
#     openssl-dev \
#     build-base \
#     python3-dev \
#     linux-headers \
#     pcre-dev \
#     gcompat


RUN python -m venv /py && \
    /py/bin/pip install --upgrade pip

# export scripts and virtualenv to path
ENV PATH="/scripts:/py/bin:$PATH"
    
RUN mkdir /barcode_identifier_api
WORKDIR /barcode_identifier_api
COPY requirements.txt /barcode_identifier_api
RUN pip install -r requirements.txt

EXPOSE 8000

COPY ./scripts /scripts

RUN chmod -R u+x /scripts && \
    addgroup --system appgroup && \
    adduser --ingroup appgroup appuser --disabled-password --gecos "" --no-create-home && \
    chown -R appuser:appgroup /scripts && \
    mkdir -p /var/www/runs && \ 
    chmod -R 764 /var/www/runs && \
    chown -R appuser:appgroup /var/www/runs && \
    # give access to runs folder to store non-served data
    mkdir -p /var/data/runs && \
    chmod -R 764 /var/data && \
    chown -R appuser:appgroup /var/data && \
    # give access to fishdb folder to store databases
    mkdir -p /var/data/fishdb && \
    chmod -R 764 /var/data/fishdb && \
    chown -R appuser:appgroup /var/data/fishdb && \
    # make sure our scripts can run
    mkdir -p /barcode_identifier_api/scripts && \
    chmod -R u+x /barcode_identifier_api/scripts && \
    chown -R appuser:appgroup /barcode_identifier_api/scripts && \
    # TODO: only allow the appuser to run the ncbi scripts
    # prepare folders to move static and media content
    mkdir -p /vol/web/static && \
    mkdir -p /vol/web/media && \
    chown -R appuser:appgroup /vol && \
    chmod -R 755 /vol

USER appuser

# FROM pyvenv_image


# WORKDIR /barcode_identifier_api

# # run celery workers as appuser
# USER appuser

# CMD [ "celery_run.sh" ] 

