#!/bin/bash

RED_COLOR_ESCAPE='\033[0;31m'
GREEN_COLOR_ESCAPE='\033[0;32m'
BLUE_COLOR_ESCAPE='\033[0;34m'
NO_COLOR_ESCAPE='\033[0m'
YELLOW_COLOR_ESCAPE='\033[1;33m'

BARREL_VERSION='0.0.3'

echo -e "Welcome to ${BLUE_COLOR_ESCAPE} Barrel v${BARREL_VERSION}"

# DB_NAME: Name of postgres database
# DB_USER: Username of postgres adminstrator
# DB_EMAIL: Email of postgres adminstrator 
# SECRET_KEY: Django secret key 
var_names=(DB_NAME DB_USER DB_PASSWORD DB_EMAIL SECRET_KEY DEBUG RABBITMQ_USERNAME RABBITMQ_PASSWORD)
var_defaults=(barcode_identifier_db admin 1234 admin@example.com SECRET_KEY_!\@\$123 0 username password123)

echo -e "${RED_COLOR_ESCAPE}WARNING: Proceeding with this setup script will overwrite any existing changes in the .env file, if one already exists from a previous installation.${NO_COLOR_ESCAPE}" 
read -p "Proceed? (y/n): " answer
if [ "$answer" == "${answer#[Yy]}" ]; then 
    echo "Exiting ..."
    exit
fi

echo "BARREL_VERSION=$BARREL_VERSION" > .env

# Prompt the user whether the server should be run locally or outside
while true;
do 
    read -r -p "Setup server to accept connections over the web? Setup for outside connections will require you to provide a valid SSL/TLS certificate. (y/n)" yn
    case $yn in 
        [Yy]* ) ENABLE_CONNECTIONS=true; break;;
        [Nn]* ) ENABLE_CONNECTIONS=false; break;;
    esac
done 

if [ $ENABLE_CONNECTIONS = true ]; then 
    echo "  ~ Barrel will setup server to accept outside connections over HTTPS. Ports 80 and 443 will be listen for incoming HTTP and HTTPS web requests."
    read -p "   Input your website name (e.g. barrel.com): " SERVER_NAME
    read -p "   Input full path to your certificate file. (e.g. /etc/ssl/certs/my_certificate.crt): " SSL_SERVER_CERT_PATH
    read -p "   Input full path to the private key file. (e.g. /etc/ssl/certs/my_certificate.key): " SSL_PRIVATE_KEY_PATH
    EXPOSED_PORT=80
    SECURE_PORT=443
    CONF_FILE="https_default.conf.tpl"
    ALLOWED_HOSTS=\*
else 
    echo "  ~ Barrel will setup server to NOT accept outside connections. Barrel will be run locally on localhost on a specified port. It will listen for HTTPS requests on another port, but requests will be error or be denied."
    SERVER_NAME=localhost 
    SSL_SERVER_CERT_PATH='./dummy.txt'
    SSL_PRIVATE_KEY_PATH='./dummy.txt'
    read -p "   Specify the port to be used (e.g. 8000): " EXPOSED_PORT 
    CONF_FILE="http_default.conf.tpl"
    ALLOWED_HOSTS=localhost
    SECURE_PORT=8001
fi

echo "CONF_FILE=$CONF_FILE" >> .env; 
echo "ENABLE_CONNECTIONS=$ENABLE_CONNECTIONS" >> .env; 
echo "SSL_SERVER_CERT_PATH=$SSL_SERVER_CERT_PATH" >> .env; 
echo "SSL_PRIVATE_KEY_PATH=$SSL_PRIVATE_KEY_PATH" >> .env; 
echo "SERVER_NAME=$SERVER_NAME" >> .env
echo "EXPOSED_PORT=$EXPOSED_PORT" >> .env 
echo "SECURE_PORT=$SECURE_PORT" >> .env 
echo "ALLOWED_HOSTS=$ALLOWED_HOSTS" >> .env 

# Supply redirect for HTTP -> HTTPS requests here, to avoid $ escape errors
echo "REDIRECT_URI=\"https://\$SERVER_NAME\\\$request_uri\"" >> .env;

for i in "${!var_names[@]}"; 
do 
    var_name=${var_names[$i]}
    read -p "   Enter value for $var_name (default: ${var_defaults[$i]}): " var_value
    # If value is empty.
    if [ -z "$var_value" ]; then 
        var_value=${var_defaults[$i]}
        echo "  ~ using default ${var_value}"
    fi
    echo "$var_name='$var_value'" >> .env
done

echo "DJANGO_SUPERUSER_USERNAME=\$DB_USER" >> .env
echo "DJANGO_SUPERUSER_EMAIL=\$DB_EMAIL" >> .env
echo "DJANGO_SUPERUSER_PASSWORD=\$DB_PASSWORD" >> .env
echo "POSTGRES_USER=\$DB_USER" >> .env
echo "POSTGRES_DB=\$DB_NAME" >> .env
echo "POSTGRES_PASSWORD=\$DB_PASSWORD" >> .env
echo "DB_HOST=barrel_db" >> .env
echo "CELERY_BROKER_URL=amqp://\$RABBITMQ_USERNAME:\$RABBITMQ_PASSWORD@barrel_rabbitmq:5672//" >> .env


echo -e "${GREEN_COLOR_ESCAPE} (server_setup.sh) Successfully created .env file in this directory. ${NO_COLOR_ESCAPE}"

