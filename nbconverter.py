import glob
import os
import os.path
import warnings
from nbconvert.exporters import markdown

warnings.simplefilter('ignore')

notebook_source_dir = '../notebooks'


def nb_to_markdown(nb_path):
    """Convert notebook to markdown"""
    exporter = markdown.MarkdownExporter()
    output, resources = exporter.from_file(open(nb_path))
    base_url = 'http://nbviewer.jupyter.org/github/unidata/notebook-gallery/tree/master/notebooks/'
    out_lines = ['[Notebook](%s)' % (base_url + os.path.basename(nb_path))]
    for line in output.split('\n'):
        out_lines.append(line)
    output = '\n'.join(out_lines)

    return output, resources


def write_nb(dest, output, resources):
    if not os.path.exists(dest):
        os.makedirs(dest)
    md_file = os.path.join(dest, resources['metadata']['basename'] + resources['output_extension'])
    name = resources['metadata']['name']
    with open(md_file, 'w') as md:
        md.write('#' + name + '\n')
        md.write(output)

    imgdir = os.path.join(dest, resources['metadata']['basename'])
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)
    basename = resources['metadata']['basename']
    for filename in resources['outputs']:
        img_file = os.path.join(imgdir, filename.replace('output_', basename + ' '))
        with open(img_file, 'wb') as img:
            img.write(resources['outputs'][filename])


def generate(app):
    for fname in glob.glob(os.path.join(app.srcdir, notebook_source_dir, '*.ipynb')):
        write_nb(os.path.join(app.srcdir, notebook_source_dir, *nb_to_markdown(fname)))
    with open(os.path.join(app.srcdir, notebook_source_dir, 'test.md'), 'w') as test:
        test.write('#Notebook Gallery\n')
        for fname in glob.glob(os.path.join(app.srcdir, notebook_source_dir, '*.ipynb')):
            filepath, filename = os.path.split(fname)
            test.write('![' + filename + '](link to picture)('
                       'http://nbviewer.jupyter.org/github/unidata/notebook-gallery/tree/master/notebooks/' +
                       filename + ')\n')
