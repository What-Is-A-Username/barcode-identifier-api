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

    location /api {
        rewrite                 ^/api/(.*)$ /$1 break;
        uwsgi_pass              ${APP_HOST}:${APP_PORT};
        include                 /etc/nginx/uwsgi_params;
        client_max_body_size    10M;
    }

    location /app {
        alias /frontend;
        index index.html;
    }
}