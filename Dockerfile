FROM python:3.8.10-alpine3.13
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1
#TODO: postgresql-client apk needed?
#TODO: Remove .tmp-deps after install
RUN apk add --no-cache postgresql-dev \ 
    gcc \
    musl-dev \
    curl-dev \
    openssl-dev \
    build-base \
    python3-dev \
    linux-headers \
    pcre-dev

RUN python -m venv /py && \
    /py/bin/pip install --upgrade pip

# export scripts and virtualenv to path
ENV PATH="/scripts:/py/bin:$PATH"
    
RUN mkdir /barcode_identifier_api
WORKDIR /barcode_identifier_api

RUN addgroup -S appgroup && \
    adduser -S appuser -g appgroup && \
    mkdir -p /var/www/runs && \ 
    chown -R appuser:appgroup /var/www/runs && \
    chmod -R 755 /var/www/runs && \
    mkdir -p /barcode_identifier_api/runs && \
    mkdir -p /barcode_identifier_api/scripts && \
    chmod -R +x /barcode_identifier_api/scripts && \
    mkdir -p /vol/web/static && \
    mkdir -p /vol/web/media && \
    chown -R appuser:appgroup /vol && \
    chmod -R 755 /vol

COPY requirements.txt /barcode_identifier_api
RUN pip install -r requirements.txt

EXPOSE 8000

COPY ./scripts /scripts

RUN chmod -R +x /scripts

USER appuser

CMD [ "run.sh" ]

