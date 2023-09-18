FROM barcode-identifier-api-pyvenv_image:latest
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1
#TODO: postgresql-client apk needed?
#TODO: Remove .tmp-deps after install

WORKDIR /barcode_identifier_api

EXPOSE 8000

USER appuser

CMD [ "run.sh" ]

