import glob
import os
import os.path
from nbconvert.exporters import markdown

SITE_DIR = 'website'
IMAGE_DIR = os.path.join(SITE_DIR, 'images')

def make_thumbnail(nb_path):
    """Convert notebook to markdown"""
    exporter = markdown.MarkdownExporter()
    output, resources = exporter.from_filename(nb_path)
    images = sorted(resources['outputs'])
    if images:
        img_file = images[-1]
        thumb_name = os.path.join(IMAGE_DIR, os.path.splitext(os.path.basename(nb_path))[0] +
                os.path.splitext(img_file)[-1])

        with open(thumb_name, 'wb') as thumb:
            thumb.write(resources['outputs'][img_file])

        return os.path.relpath(thumb_name, SITE_DIR)
    else:
        return 'images/placeholder.png'


if __name__ == '__main__':
    with open(os.path.join(SITE_DIR, 'index.md'), 'w') as index:
        index.write('---\ntitle: Unidata\'s Notebook Gallery\n---\n# Notebook Gallery\n\n')
        base_url = 'http://nbviewer.jupyter.org/github/unidata/notebook-gallery/blob/master/'
        index.write('<div id="gallery">\n')
        for fname in sorted(glob.glob(os.path.join('notebooks', '*.ipynb'))):
            print('Converting {0} -> '.format(fname), end='')
            img_file = make_thumbnail(fname)
            print(img_file)
            index.write('<div><a href="{url}{nb_file}"><img src="{img}"></a></div>\n'.format(url=base_url,
                nb_file=fname, img=img_file))
        index.write('</div>\n')
