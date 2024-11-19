FROM python:3.8.10
LABEL Author="barrel_authors"

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
    libffi-dev \
    libmagic1

RUN python -m venv /py && \
    /py/bin/pip install --upgrade pip

RUN mkdir /barrel
WORKDIR /barrel

# export scripts and virtualenv to path
ENV PATH="/barrel/scripts:/py/bin:$PATH"

COPY ./requirements.txt /barrel/
RUN pip install -r requirements.txt

COPY ./barcode_identifier_api /barrel/barcode_identifier_api
COPY ./barcode_blastn /barrel/barcode_blastn
COPY ./manage.py /barrel/ 
COPY ./docker.env /barrel/
COPY ./scripts /barrel/scripts
COPY ./compose /barrel/compose
COPY ./ncbi-blast-2.12.0+ /barrel/ncbi-blast-2.12.0+

RUN chmod -R u+x /barrel/scripts && \
    addgroup --system appgroup && \
    adduser --ingroup appgroup appuser --disabled-password --gecos "" --no-create-home && \
    chown -R appuser:appgroup /barrel/scripts && \
    chown -R appuser:appgroup /barrel/barcode_blastn/migrations && \
    touch /barrel/barcode_blastn/migrations/__init__.py && \
    # Ensure that the soon-to-be named volumes have the correct permissions for other containers
    mkdir -p /vol/static && \
    mkdir -p /var/data && \
    chmod 755 /vol/static && \
    chmod 755 /var/data && \
    chown -R appuser:appgroup /vol/static && \ 
    chown -R appuser:appgroup /var/data
    # mkdir -p /var/www/runs && \ 
    # chmod -R 764 /var/www/runs && \
    # chown -R appuser:appgroup /var/www/runs && \
    # give access to runs folder to store non-served data
    # mkdir -p /var/data/runs && \
    # chmod -R 764 /var/data && \
    # chown -R appuser:appgroup /var/data && \
    # # give access to fishdb folder to store databases
    # mkdir -p /var/data/library && \
    # chmod -R 764 /var/data/library && \
    # chown -R appuser:appgroup /var/data/library
    # make sure our scripts can run
    # mkdir -p /barrel/scripts && \
    # chmod -R u+x /barrel/scripts && \
    # chown -R appuser:appgroup /barrel/scripts
    # TODO: only allow the appuser to run the ncbi scripts
    # prepare folders to move static and media content
    # mkdir -p /vol/web/static && \
    # mkdir -p /vol/web/media && \
    # chown -R appuser:appgroup /vol && \
    # chmod -R 755 /vol

EXPOSE 8000

