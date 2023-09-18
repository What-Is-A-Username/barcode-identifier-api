import magic
from django.template.defaultfilters import filesizeformat
from django.forms import ValidationError

class QueryFileValidator(object):
    default_error_messages = {
        'max_size': ("The uploaded file cannot be greater than %(max_size)s."
                    " The present file size is %(size)s."),
        'content_type': "The uploaded file type %(content_type)s cannot be used for the query file.",
    }

    def __init__(self, max_size: int = 3145728):
        self.max_size = max_size # 3 MB max file size
        self.content_types = ['text/plain']

    def __call__(self, data):
        if self.max_size is not None and data.size > self.max_size:
            raise ValidationError(self.default_error_messages['max_size'], 
                'max_size', 
                {
                    'max_size': filesizeformat(self.max_size), 
                    'size': filesizeformat(data.size),
                })
        if self.content_types:
            content_type = magic.from_buffer(data.read(), mime=True)
            data.seek(0)

            if content_type not in self.content_types:
                raise ValidationError(self.default_error_messages['content_type'], 
                'content_type', 
                { 'content_type': content_type })

            return data
        else:
            return data

    def __eq__(self, other):
        return (
            isinstance(other, QueryFileValidator) and
            self.max_size == other.max_size and
            self.content_types == other.content_types
        )

