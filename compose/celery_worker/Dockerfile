FROM barcode_identifier_api_pyvenv_image:latest
LABEL Author="new_author here"

ENV PYTHONBUFFERED 1   

WORKDIR /barcode_identifier_api

# RUN addgroup -S celerygroup && \
#     adduser -S celeryuser -g celerygroup
    # give access to runs folder to store static content to be served
# RUN ls -l /var/ && \
#     mkdir -p /var/www/runs && \ 
#     chmod -R 764 /var/www/runs && \
#     chown -R appuser:appgroup /var/www/runs && \
#     # give access to runs folder to store non-served data
#     mkdir -p /var/data/runs && \
#     chmod -R 764 /var/data/runs && \
#     chown -R appuser:appgroup /var/data/runs && \
#     # give access to fishdb folder to store databases
#     mkdir -p /var/data/fishdb && \
#     chmod -R 764 /var/data/fishdb && \
#     chown -R appuser:appgroup /var/data/fishdb && \
#     # make sure our scripts can run
#     mkdir -p /barcode_identifier_api/scripts && \
#     chmod -R u+x /barcode_identifier_api/scripts && \
#     chown -R appuser:appgroup /barcode_identifier_api/scripts && \
#     # TODO: only allow the appuser to run the ncbi scripts
#     # prepare folders to move static and media content
#     mkdir -p /vol/web/static && \
#     mkdir -p /vol/web/media && \
#     chown -R appuser:appgroup /vol && \
#     chmod -R 755 /vol

# run celery workers as appuser
USER appuser

# CMD [ "celery_run.sh" ] 

