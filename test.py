# uwsgi test, Nov 27th 2022, according to https://uwsgi-docs.readthedocs.io/en/latest/tutorials/Django_and_nginx.html#basic-test

# test.py
def application(env, start_response):
    start_response('200 OK', [('Content-Type','text/html')])
    return [b"Hello World"] # python3
    #return ["Hello World"] # python