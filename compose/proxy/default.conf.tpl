server {
    listen ${LISTEN_PORT};
    server_name ${SERVER_NAME};

    location /static/runs/ {
        deny all;
        alias /vol/static/runs/;
        # only allow .clustal_num, .ph, and .txt files
        location ~ \.(clustal_num|ph|txt)$ {
            add_header X-Robots-Tag "noindex, follow" always;
	        allow all;
        }
    }

    location /static {
        alias /vol/static;
        add_header X-Robots-Tag "noindex, follow" always;
    }

    location /api {
        uwsgi_pass              ${APP_HOST}:${APP_PORT};
        include                 /etc/nginx/uwsgi_params;
        client_max_body_size    10M;
        add_header X-Robots-Tag "noindex, follow" always;
    }

    location /app/assets {
        alias /frontend/assets;
        add_header X-Robots-Tag "noindex, follow" always;
    }

    location /app {
        alias /frontend;
        location ~* \.(jpeg|jpg|png|txt|json|svg|ico)$ {
            expires 1h;
            add_header X-Robots-Tag "noindex, follow" always;
        }
        location /app {
            try_files /index.html =404;
            add_header X-Robots-Tag "noindex, follow" always;
        }
    }
}