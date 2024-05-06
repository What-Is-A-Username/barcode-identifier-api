FROM nginxinc/nginx-unprivileged:1-alpine
#TODO: Edit maintainer to email or website link
LABEL maintainer="barcode_identifier_dev"

COPY ./uwsgi_params /etc/nginx/uwsgi_params 
COPY --chown=nginx:nginx ./run.sh /run.sh

# Copy SSL certificates for SSL
ARG SSL_SERVER_CERT_PATH
ENV SSL_SERVER_CERT_PATH $SSL_SERVER_CERT_PATH
ARG SSL_PRIVATE_KEY_PATH
ENV SSL_PRIVATE_KEY_PATH $SSL_PRIVATE_KEY_PATH
COPY --chown=root:nginx $SSL_SERVER_CERT_PATH /etc/ssl/certs/barrel.crt
COPY --chown=root:nginx $SSL_PRIVATE_KEY_PATH /etc/ssl/private/barrel.key

# COPY --chown=root:nginx ./$SSL_CONF_FILE /etc/ssl/ssl.conf

# Name of service running Django app
ENV APP_HOST=barrel
# Barrel is running on port 9000
ENV APP_PORT=9000

# Based on the setup options, use appropriate server config
ARG CONF_FILE
COPY ./${CONF_FILE} /etc/nginx/default.conf.tpl 

USER root

# make folder for static 
RUN rm -f /etc/nginx/conf.d/default.conf && \
    # make empty file 
    touch /etc/nginx/conf.d/default.conf && \
    # change ownership to nginx user 
    chown nginx:nginx /etc/nginx/conf.d/default.conf && \
    # allow script to execute
    chmod 700 /run.sh && \
    # restrict permissions to protect certificate and key
    # only allow the nginx user to read
    chmod 644 /etc/ssl/certs/barrel.crt && \
    chmod 640 /etc/ssl/private/barrel.key && \
    chown root:nginx -R /etc/ssl/certs/barrel.crt && \
    chown root:nginx -R /etc/ssl/private/barrel.key 

# switch to nginx user
USER nginx 

CMD ["./run.sh"]