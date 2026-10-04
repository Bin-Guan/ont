#!/bin/bash
#$1 version on GitHub, the git tag command does not match with GitHub.com??
#$2 message, quotation mark for sentence
echo $1 "$2"
set -euo pipefail
git add --all
git commit -m "$2"
git tag --force -a "$1" -m "$2"
git log -1 > /data/OGL/resources/ont.git.log
printf 'Tag: %s\n' "$1" >> /data/OGL/resources/ont.git.log
branch=$(git branch --show-current)
git push --set-upstream origin "$branch"
git push --force origin "$1"

# git add --all
# git commit -a -m "$2"
# git tag --force -a $1 -m "$2"
# git log | head -n 5 > /data/OGL/resources/ont.git.log
# git tag | tail -n 1 >> /data/OGL/resources/ont.git.log
# branch=$(git branch --show-current)
# git push --set-upstream origin "$branch"
# git push origin "$1"