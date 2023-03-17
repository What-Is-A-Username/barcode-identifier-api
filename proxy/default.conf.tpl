server {
    listen ${LISTEN_PORT};

    location /static/runs/ {
        deny all;
        alias /var/www/runs/;
        # only allow .clustal_num, .ph, and .txt files
        location ~ \.(clustal_num|ph|txt)$ {
            add_header X-Robots-Tag "noindex, follow" always;
	        allow all;
        }
    }

    location /static {
        alias /vol/static;
    }

    location / {
        {* where APP_HOST is container containing app, APP_PORT is the port on it *}
        uwsgi_pass              ${APP_HOST}:${APP_PORT};
        include                 /etc/nginx/uwsgi_params;
        {* max size of request is 10 MB, to limit file upload*}
        client_max_body_size    10M;
    }
}