#!/bin/bash
set -e # exit with nonzero exit code if anything fails

# Decrypt and activate the deploy key
echo Setting up access...
openssl aes-256-cbc -K $encrypted_0e2eca647d99_key -iv $encrypted_0e2eca647d99_iv -in deploy/deploy_key.enc -out deploy/deploy_key -d
chmod 600 deploy/deploy_key
eval `ssh-agent -s`
ssh-add deploy/deploy_key

# inside this git repo we'll pretend to be a new user

echo Copying...
cp -R ${TRAVIS_BUILD_DIR}/website/build/html $HOME/gh_pages
touch .nojekyll

echo Staging...
cd $HOME/gh_pages
git init
git config user.name "Travis CI"
git config user.email "travis@nobody.org"
git add -A .
git commit -m "Deploy to GitHub Pages"

# Push up to gh-pages
echo Pushing...
git push --force git@github.com:${TRAVIS_REPO_SLUG}.git master:gh-pages
