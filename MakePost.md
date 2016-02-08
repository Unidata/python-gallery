# Turning a notebook into a post

To turn a notebook into a post, from the root of the repository (where the
specified python file lives) run the following:

```sh
jupyter nbconvert --config blog-post.py --to markdown notebooks/NOTEBOOK_FILE.ipynb
```

This will generate a `.md` file with the contents that can be copied/pasted
into the blog post editor. On osx, it is even easier:

```sh
jupyter nbconvert --config blog-post.py --to markdown --stdout notebooks/NOTEBOOK_FILE.ipynb | pbcopy
```

This will convert the notebook and put the content into the paste buffer.
