FROM barcode-identifier-api-barrel_venv_image:latest
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1
#TODO: postgresql-client apk needed?
#TODO: Remove .tmp-deps after install

WORKDIR /barrel

EXPOSE 9000

USER appuser

CMD [ "./scripts/run.sh" ]

