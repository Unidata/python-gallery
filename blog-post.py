try:
    from urllib.parse import quote  # Py 3
except ImportError:
    from urllib2 import quote  # Py 2
import os
import sys
from datetime import datetime

today = datetime.utcnow()

c = get_config()
c.NbConvertApp.export_format = 'markdown'
c.MarkdownExporter.template_file = 'blog'
c.ExtractOutputPreprocessor.output_filename_template = today.strftime('%Y%m%d') + '_{unique_key}_{cell_index}_{index}{extension}'

def path2url(path):
    """Turn a file path into a URL"""
    parts = path.split(os.path.sep)
    img_dir = '/blog_content/images/%s/' % today.strftime('%Y')
    return img_dir + '/'.join(quote(part) for part in parts[1:])

c.MarkdownExporter.filters = {'path2url': path2url}
