FROM barrel-barrel_venv_image:latest
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1
#TODO: postgresql-client apk needed?
#TODO: Remove .tmp-deps after install

WORKDIR /barrel

EXPOSE 9000
COPY ./scripts/run.sh /barrel/scripts/run.sh
RUN ls -lha /barrel/scripts/ && \
    chmod -R u+x /barrel/scripts/run.sh && \
    chown -R appuser:appgroup /barrel/scripts/run.sh

USER appuser

# CMD [ "./wow/run.sh" ]

