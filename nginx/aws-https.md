# Setting up HTTPS using openssl cookbook (obselete)

I am following openssl book here: https://www.feistyduck.com/library/openssl-cookbook/online/openssl-command-line/key-generation.html. 

## Key generation
Generate the `fd.key` and set the PEM pass phrase when prompted.
```
openssl genpkey -out fd.key -algorithm RSA -pkeyopt rsa_keygen_bits:2048 -aes-128-cbc
```

The key generation page of the book also provides other generation methods to generate fd.key.

## Create certificate signing request
```
openssl req -new -key fd.key -out fd.csr
```

# Setting up HTTPS using uWSGI docs 

Note; *The following procedure generates a self-signed certificate to enable HTTPS. Most browsers will not trust self-signed certificate, and may prevent you from accessing the page*. **Use only for development, not production.**

The entry in the docs is at https://uwsgi.readthedocs.io/en/latest/HTTPS.html

```
cd ~
mkdir ~/cert
cd ~/cert 
openssl genrsa -out foobar.key 2048
openssl req -new -key foobar.key -out foobar.csr
openssl x509 -req -days 365 -in foobar.csr -signkey foobar.key -out foobar.crt
sudo cp ~/cert/** /etc/nginx/ssl/
```

Open `barcode_identifier_api_nginx.conf`,
```
sudo nano /etc/nginx/sites-available/barcode_identifier_api_nginx.conf
```
and update it to look something like this:
```
# copy to /etc/nginx/sites-available/barcode_identifier_api_nginx.conf
# the upstream component nginx needs to connect to
upstream django {
    server unix:///home/ubuntu/barcode_identifier_api/barcode_identifier_api.sock;
    # server 127.0.0.1:8001; # for a web port socket (we'll use this first)
}

server {
    listen      80;
    server_name     3.17.163.210;
    rewrite ^/(.*)  https://3.17.163.210/$1 permanent;
}

# configuration of the server
server {
    # the port your site will be served on
    listen      443 ssl;
    # the domain name it will serve for
    server_name 3.17.163.210; # substitute your machine's IP address or FQDN
    # server_name placeholder.com;
    charset     utf-8;

    # ssl certificate for https
    ssl_certificate     /etc/nginx/ssl/foobar.crt;
    ssl_certificate_key     /etc/nginx/ssl/foobar.key;

    # max upload size
    client_max_body_size 75M;   # adjust to taste

    # Django media
    location /media  {
        alias /home/ubuntu/barcode_identifier_api/media;
    }

    location /static {
        alias /home/ubuntu/barcode_identifier_api/static;
    }

    # Finally, send all non-media requests to the Django server.
    location / {
        uwsgi_pass  django;
        include     /home/ubuntu/barcode_identifier_api/uwsgi_params;
    }
}
```

# Setting up HTTPS with Let's Encrypt and Certbot (incomplete)

The following information is the first few steps for the process of obtaining and installing an SSL certificate for the domain through "Lets Encrypt" and "Certbot". **SSL certificate only obtainable for a non amazonaws.com domain.**

```
sudo snap install --classic certbot
sudo ln -s /snap/bin/certbot /usr/bin/certbot
```
