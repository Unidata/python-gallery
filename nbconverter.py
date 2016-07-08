import glob
import os
import os.path
import warnings
from nbconvert.exporters import markdown

warnings.simplefilter('ignore')

notebook_source_dir = 'notebooks'


def nb_to_markdown(nb_path):
    """Convert notebook to markdown"""
    exporter = markdown.MarkdownExporter()
    output, resources = exporter.from_file(open(nb_path))
    basename = os.path.splitext(os.path.basename(nb_path))[0]
    resources['metadata']['basename'] = basename
    base_url = 'http://nbviewer.jupyter.org/github/unidata/notebook-gallery/tree/master/notebooks/'
    out_lines = ['[Notebook](%s)' % (base_url + os.path.basename(nb_path))]
    for line in output.split('\n'):
        out_lines.append(line)
    output = '\n'.join(out_lines)

    return output, resources


def write_nb(dest, output, resources):
    if not os.path.exists(dest):
        os.makedirs(dest)
    md_file = os.path.join('website', resources['metadata']['basename'] + resources['output_extension'])
    name = resources['metadata']['name']
    with open(md_file, 'w') as md:
        md.write('---\ntitle:\n--- \n#' + name + '\n')
        md.write(output)

    imgdir = os.path.join('website', resources['metadata']['basename'])
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)
    basename = resources['metadata']['basename']
    img_file = ""
    for filename in resources['outputs']:
        img_file = os.path.join(imgdir, filename.replace('output_', basename))
        with open(img_file, 'wb') as img:
            img.write(resources['outputs'][filename])
    return img_file


if __name__ == '__main__':
    with open(os.path.join('website', 'index.md'), 'w') as test:
        test.write('---\ntitle: Unidata\'s Notebook Gallery\n--- \n#Notebook Gallery\n')
        for fname in glob.glob(os.path.join(notebook_source_dir, '*.ipynb')):
            img_file = write_nb(notebook_source_dir, *nb_to_markdown(fname))
            filename = os.path.split(fname)[1]
            if not img_file:
                img_file = 'website/images/placeholder.png'
            test.write('<a href="notebooks/' + filename + '"><img alt="' + filename.replace('_', ' ').replace('.ipynb', '') +
                       '"src="' + img_file + '" height="300" width="375"></a>\n')
