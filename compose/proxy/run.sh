#!/bin/sh

# 
set -e

# move config into default.conf while substituting in ENV variables
envsubst < /etc/nginx/default.conf.tpl > /etc/nginx/conf.d/default.conf
# start in nginx in the foreground, not the background, so debug is sent to log
nginx -g 'daemon off;'